function results = MANSS(netlistFilePath, config)
%AnalyzeSwitchingConverter Performs PSS analysis and ZVS deadtime
%optimization. 

%    DO NOT MODIFY - RUNNER SCRIPT (run_MANSS.m) CONTROLS

%    SWITCH PAIRING RECOGNITION ADDED
%    CUSTOM PARASITICS ADDED
%
%   results = AnalyzeSwitchingConverter(netlistFilePath, config)
%
%   This function takes a SPICE-like netlist file and a configuration
%   structure to perform a comprehensive analysis of a switching power
%   converter. It automatically:
%   1. Parses the netlist (R, L, C, V, I, E, G, H, F, K, M).
%   2. Determines valid switching states, excluding shoot-through.
%   3. Generates state-space matrices (A, B, C, D) for each valid state using
%      Modified Nodal Analysis (MNA).
%      Handles mutual inductance (K). Detects MNA degeneracy and errors out.
%   4. Includes optional RC switch parasitic modeling (Ron/Roff || (Rs + Cs)).
%   5. Calculates PSS initial conditions using Enhanced SVA (ESVA), based on Sadri (2022).
%      Handles potentially singular (I - Phi_composite) using pinv with warning.
%   6. Optimizes deadtimes for ZVS using Newton-Raphson, based on logic from
%      Ellouz_ESVA_Analysis_NetlistImport_Update9_4.txt
%   7. Generates steady-state waveforms for state variables and specified outputs.
%   8. Returns a structure containing the analysis results.
%
%   Inputs:
%       netlistFilePath - String: Path to the netlist file.
%       config          - Struct: Contains configuration parameters.
%                         Required fields: Ron, Roff, fs, D.
%                         See runner script for all options.
%
%   Outputs:
%       results         - Struct: Contains analysis results (x0_ss, optimal_dt, waveforms, etc.).
%
% 
%
%   Version: 2.0 (2025-04-17) - Full implementation based on detailed plan.

% =========================================================================
%                       Input Validation & Setup
% =========================================================================
    if nargin < 2
        error('AnalyzeSwitchingConverter:NotEnoughInputs', ...
              'Requires netlistFilePath (string) and config (struct) as input.');
    end
    if ~ischar(netlistFilePath) && ~isstring(netlistFilePath)
        error('AnalyzeSwitchingConverter:InvalidNetlistPathType', 'netlistFilePath must be a string.');
    end
     if ~isfile(netlistFilePath)
        error('AnalyzeSwitchingConverter:NetlistNotFound', 'Netlist file not found: %s', netlistFilePath);
    end
    if ~isstruct(config)
         error('AnalyzeSwitchingConverter:InvalidConfigType', 'config must be a struct.');
    end

    % Start overall timer
    timerVal = tic;
    fprintf('\n========== Starting AnalyzeSwitchingConverter ==========\n');

    % Use a main try-catch block for robust error handling
    try
% =========================================================================
%                         Phase 1: Initialization & Input Processing
% =========================================================================
fprintf('\n--- Phase 1: Initialization & Input Processing ---\n');

% --- Validate and Merge Configuration ---
 config = validateAndMergeConfig(config); % Assuming config is already validated by runner


% --- Parse Netlist ---
fprintf('Parsing netlist file: %s\n', netlistFilePath);
% Added parsed_input_vals output from parseNetlist
[components, switches, node_map_initial, coupling_map, config, u_vars, parsed_input_vals] = parseNetlist(netlistFilePath, config); % Nested function
% --- Robustly define num_switches based on the 'switches' variable from parseNetlist ---
    % The 'switches' variable below is the direct output from the parseNetlist call above.
    
    fprintf('DEBUG: Immediately after parseNetlist call. Checking ''switches'' variable.\n');
    if ~exist('switches', 'var')
        % This case should ideally not be reached if 'switches' is an output argument of parseNetlist
        % and parseNetlist itself doesn't error out before assigning outputs.
        fprintf(2, 'CRITICAL DEBUG ERROR: ''switches'' variable DOES NOT EXIST in the workspace immediately after the call to parseNetlist.\n');
        fprintf(2, 'This implies a fundamental issue with parseNetlist not returning it or a workspace problem.\n');
        % Failsafe: Initialize 'switches' as an empty struct array with expected fields
        switches = struct('name', {}, 'type', {}, 'nodes', {}, 'ctrl_node', {}, 'index', {}); 
        num_switches = 0; 
        fprintf(2, 'Failsafe: ''switches'' has been re-initialized as an empty struct, num_switches set to 0.\n');
    else
        % 'switches' variable exists. Proceed to check its type and determine num_switches.
        fprintf('DEBUG: ''switches'' variable EXISTS after parseNetlist call. Class: %s, Size: %s\n', class(switches), mat2str(size(switches)));
        if ~isempty(switches)
            % whos switches; % Uncomment for very detailed info if needed
            fprintf('DEBUG: ''switches'' variable is not empty.\n');
        else
            fprintf('DEBUG: ''switches'' variable is empty.\n');
        end
        
        try
            if isstruct(switches)
                % If switches = struct('name', {}) (an empty struct with fields), length is 0.
                % If switches = struct([]) (an empty struct with no fields), length is 1, which is incorrect for num_switches.
                if isempty(switches) && ~isfield(switches,'name') % Specifically handles struct([])
                    num_switches = 0;
                    fprintf('DEBUG: ''switches'' is struct([]), so num_switches set to 0.\n');
                else % Handles non-empty struct array AND empty struct array with fields (e.g. struct('name',{}))
                    num_switches = length(switches);
                    fprintf('DEBUG: num_switches defined as %d using length(switches) for struct.\n', num_switches);
                end
            elseif iscell(switches) && isempty(switches)
                num_switches = 0; % Handle if parseNetlist returns an empty cell for no switches
                fprintf('DEBUG: ''switches'' is an empty cell, num_switches set to 0.\n');
            else 
                fprintf(2, 'WARNING: ''switches'' output from parseNetlist is of unexpected type: %s. Assuming 0 switches.\n', class(switches));
                switches = struct('name', {}, 'type', {}, 'nodes', {}, 'ctrl_node', {}, 'index', {}); % Failsafe reset
                num_switches = 0; % Failsafe
            end
        catch ME_len
            fprintf(2, 'DEBUG ERROR: Could not execute length(switches) or other operations on ''switches''. Error: %s\n', ME_len.message);
            num_switches = -1; % Indicate an error state for num_switches
            fprintf(2, 'DEBUG: num_switches set to -1 (error indicator).\n');
        end
    end
    
    % Final check to ensure num_switches is a scalar non-negative integer
    if ~isscalar(num_switches) || ~isnumeric(num_switches) || num_switches < 0 || mod(num_switches,1) ~= 0
        warning('AnalyzeSwitchingConverter:InvalidNumSwitchesValue', ...
                'Value of num_switches (%g) is invalid after processing. Forcing to 0.', num_switches);
        num_switches = 0;
    end
    % --- End definition of num_switches ---
num_inputs = length(u_vars); % Number of independent V/I sources
fprintf('Parsing complete. Found %d components, %d nodes (incl. gnd), %d switches, %d inputs.\n', ...
        length(components), length(keys(node_map_initial)), num_switches, num_inputs);

% --- Validate Switch Parasitics Structure Content (after switches are known) ---
if config.include_parasitics && num_switches > 0 && ~config.stop_after_phase1
    fprintf('Validating individual switch parasitics configuration...\n');
    if ~isfield(config, 'switch_parasitics') || ~isstruct(config.switch_parasitics)
        % Check if it exists and is a struct (basic check)
         error('AnalyzeSwitchingConverter:MissingSwitchParasiticsStruct', 'config.include_parasitics is true, but config.switch_parasitics is missing or not a struct.');
    end
    parasitic_struct = config.switch_parasitics;
    missing_switches = {};
    invalid_switches = {};
    all_switch_names = {switches.name}; % Get names of all identified switches

    for k_sw_val = 1:num_switches
        sw_name_val = all_switch_names{k_sw_val};
        if ~isfield(parasitic_struct, sw_name_val)
            missing_switches{end+1} = sw_name_val;
        else
            sw_params = parasitic_struct.(sw_name_val);
            % Check if the field contains a struct with valid Cs and Rs
            if ~isstruct(sw_params) || ~isfield(sw_params, 'Cs') || ~isfield(sw_params, 'Rs') || ...
               ~isscalar(sw_params.Cs) || ~isnumeric(sw_params.Cs) || isnan(sw_params.Cs) || sw_params.Cs <= 0 || ...
               ~isscalar(sw_params.Rs) || ~isnumeric(sw_params.Rs) || isnan(sw_params.Rs) || sw_params.Rs < 0
                 % Capture potentially problematic values for error message
                 cs_val_err = NaN; rs_val_err = NaN;
                 if isstruct(sw_params) && isfield(sw_params, 'Cs') && isnumeric(sw_params.Cs); cs_val_err = sw_params.Cs(1); end
                 if isstruct(sw_params) && isfield(sw_params, 'Rs') && isnumeric(sw_params.Rs); rs_val_err = sw_params.Rs(1); end
                 invalid_switches{end+1} = sprintf('%s (Invalid Cs/Rs: Cs=%g, Rs=%g)', sw_name_val, cs_val_err, rs_val_err);
            end
        end
    end
    if ~isempty(missing_switches)
        error('AnalyzeSwitchingConverter:MissingSwitchParasitics', 'Parasitic values missing for switch(es): %s in config.switch_parasitics. Please define Cs (>0) and Rs (>=0) for all identified switches when include_parasitics is true.', strjoin(missing_switches, ', '));
    end
    if ~isempty(invalid_switches)
        error('AnalyzeSwitchingConverter:InvalidSwitchParasitics', 'Invalid parasitic values (Cs must be > 0, Rs must be >= 0) for switch(es): %s.', strjoin(invalid_switches, '; '));
    end
    fprintf('  Individual switch parasitics configuration validated for %d switches.\n', num_switches);
elseif config.include_parasitics && num_switches == 0 && ~config.stop_after_phase1
    fprintf('  Parasitics included, but no switches found in netlist. No parasitic values needed.\n');
elseif ~config.include_parasitics && ~config.stop_after_phase1
    fprintf('  Parasitics not included (config.include_parasitics is false).\n');
end

% --- Generate Valid Switching States ---
fprintf('Generating valid switching states...\n');
[valid_states, num_valid_states] = generateValidStates(switches, config); % Nested function
if isempty(valid_states) && num_switches > 0
     error('AnalyzeSwitchingConverter:NoValidStates', 'No valid switching states remain after filtering shoot-through.');
elseif isempty(valid_states) && num_switches == 0
    fprintf('  No switches found. Circuit has 1 static operating state.\n');
else
    fprintf('  Found %d valid switching states.\n', num_valid_states);
    disp('Valid States Matrix (rows are states, cols are switches):');
    disp(valid_states);
end

% --- Validate State Sequence ---
config.state_sequence = validateStateSequence(config.state_sequence, num_valid_states); % Nested function
fprintf('  Using state sequence: %s\n', mat2str(config.state_sequence));


% --- Check if stopping after Phase 1 ---
if config.stop_after_phase1
    fprintf('\n--- Stopping after Phase 1 as requested by config.stop_after_phase1 ---\n');
    % Package minimal results needed for user to proceed
    results = struct();
    results.config = config; % The minimal config used
    results.valid_states = valid_states;
    results.num_valid_states = num_valid_states;
    results.switches = switches;
    results.components = components; % Original parsed components
    results.node_map_initial = node_map_initial;
    results.u_vars = u_vars;
    results.parsed_input_values = parsed_input_vals;
    results.message = 'Analysis stopped after Phase 1. Check console for Valid States Matrix.';
    return; % Exit the function
end
fprintf('DEBUG MANSS: Received config.user_complementary_switch_pairs:\n');

% =========================================================================
%                         Phase 2: Generating State-Space Models (A, B, C, D)
% =========================================================================
fprintf('\n--- Phase 2: Generating State-Space Models (A, B, C, D) ---\n');

% Initialize storage
A_matrices = cell(1, num_valid_states); B_matrices = cell(1, num_valid_states);
C_matrices = cell(1, num_valid_states); D_matrices = cell(1, num_valid_states);
state_vars_list = cell(1, num_valid_states); component_lists = cell(1, num_valid_states);
state_dim = -1; 
first_state_vars = {}; consistent_state_vars = true;
node_map = node_map_initial; 

if config.include_parasitics && num_switches > 0
    fprintf('Updating node map for parasitic components...\n');
    node_map = updateNodeMapForParasitics(node_map, switches); 
    fprintf('  Node map updated. Total nodes: %d\n', length(keys(node_map)));
end
node_map_final = node_map; 

% --- Determine which states to process for model generation ---
indices_to_process = []; % Initialize as empty

% Priority 1: User-defined state index sequence (config.user_defined_state_index_sequence)
if isfield(config, 'user_defined_state_index_sequence') && ~isempty(config.user_defined_state_index_sequence)
    fprintf('  Selective state modeling enabled based on config.user_defined_state_index_sequence.\n');
    user_indices = config.user_defined_state_index_sequence;
    
    % Validate that all user-defined indices are valid row numbers for 'valid_states'
    if any(user_indices < 1) || any(user_indices > num_valid_states)
        error('AnalyzeSwitchingConverter:InvalidUserDefinedStateIndex', ...
              'One or more indices in config.user_defined_state_index_sequence are out of bounds for the generated valid_states matrix (1 to %d).', num_valid_states);
    end
    indices_to_process = unique(user_indices(:)'); % Ensure row vector and unique
    fprintf('  Will attempt to generate models for %d unique states from config.user_defined_state_index_sequence.\n', length(indices_to_process));

% Priority 2: User-defined operational state vectors (config.user_operational_state_vectors)
elseif isfield(config, 'user_operational_state_vectors') && ~isempty(config.user_operational_state_vectors)
    fprintf('  Selective state modeling enabled based on config.user_operational_state_vectors.\n');
    user_states_to_model = config.user_operational_state_vectors;
    num_user_states = length(user_states_to_model);
    temp_indices_to_process = zeros(1, num_user_states);
    num_found_user_states = 0;

    if num_switches == 0 && num_user_states > 0 % Should not happen if num_switches also checked for user_operational_state_vectors
        error('AnalyzeSwitchingConverter:UserStatesWithNoSwitches', ...
              'config.user_operational_state_vectors provided, but no switches were identified in the netlist.');
    end
    
    if num_switches > 0 % Only proceed if there are switches to match vectors against
        for k_user_state = 1:num_user_states
            current_user_vector = user_states_to_model{k_user_state};
            if length(current_user_vector) ~= num_switches
                error('AnalyzeSwitchingConverter:UserOpStateVectorLengthMismatch', ...
                      'User-defined operational state vector %d (length %d) does not match number of switches (%d).', ...
                      k_user_state, length(current_user_vector), num_switches);
            end
            [is_member, loc_in_valid_states] = ismember(current_user_vector, valid_states, 'rows');
            if is_member
                num_found_user_states = num_found_user_states + 1;
                temp_indices_to_process(num_found_user_states) = loc_in_valid_states;
            else
                error('AnalyzeSwitchingConverter:UserOpStateNotFoundInValidStates', ...
                      'User-defined operational state vector [%s] (entry %d in config.user_operational_state_vectors) was not found in the generated valid_states matrix.', ...
                      num2str(current_user_vector), k_user_state);
            end
        end
        if num_found_user_states > 0
            indices_to_process = unique(temp_indices_to_process(1:num_found_user_states));
            fprintf('  Will attempt to generate models for %d unique user-specified operational states (from vectors).\n', length(indices_to_process));
        else
             % This path means user_operational_state_vectors was not empty, but none matched.
             % The error inside the loop would have already triggered for the first non-match.
             % If it was empty to begin with, this 'else' block isn't hit.
            fprintf('  No valid user-specified operational states (from vectors) found. Defaulting to all %d valid states.\n', num_valid_states);
            indices_to_process = 1:num_valid_states;
        end
    else % No switches, but user_operational_state_vectors might have been non-empty (though it shouldn't be logical)
        if ~config.suppress_console_output && ~isempty(user_states_to_model)
            warning('AnalyzeSwitchingConverter:UserOpStatesWithNoSwitches', 'config.user_operational_state_vectors was provided, but no switches exist. No models will be generated based on these vectors.');
        end
        indices_to_process = 1:num_valid_states; % Which will be 1:0 if num_valid_states is 0 (no switches)
    end
else
    % Default: Process all valid states if neither user_defined_state_index_sequence nor user_operational_state_vectors is provided
    if ~config.suppress_console_output
        fprintf('  Selective state modeling not enabled by user. Processing all %d valid states (if any).\n', num_valid_states);
    end
    indices_to_process = 1:num_valid_states; % Process all if num_valid_states > 0, or 1:0 if 0.
end

if isempty(indices_to_process) && num_valid_states > 0 && num_switches > 0 
    % This implies user provided sequences/vectors, but they resulted in an empty set of valid indices to process.
    warning('AnalyzeSwitchingConverter:NoStatesToProcessAfterFilter', 'After applying user filters, no states are selected for model generation, but valid states exist. Check configuration. No models will be generated.');
    % There is no need to error out at this point here, subsequent functions (like enhancedSVA) will fail if they need models that weren't generated.
elseif isempty(indices_to_process) && num_switches == 0 && num_valid_states == 1
    % This is the case of no switches, one static state. indices_to_process might be 1:0.
    % We need to ensure the single static state (index 1) is processed.
    if num_valid_states == 1 
        indices_to_process = 1;
        if ~config.suppress_console_output
             fprintf('  No switches found. Will process the single static state (index 1).\n');
        end
    end
end


% --- Loop through each state TO BE PROCESSED ---
for k_proc_loop = 1:length(indices_to_process)
    state_loop_idx = indices_to_process(k_proc_loop); 

    current_state_vector = []; 
    if num_switches > 0
        if state_loop_idx > size(valid_states,1) || state_loop_idx < 1
            error('InternalError:InvalidStateLoopIdx', 'state_loop_idx %d is out of bounds for valid_states matrix during model generation loop.', state_loop_idx);
        end
        current_state_vector = valid_states(state_loop_idx, :);
        fprintf('Generating model for State Index %d (Vector: [%s])\n', state_loop_idx, num2str(current_state_vector));
    elseif state_loop_idx == 1 && num_switches == 0 % Only one static state if no switches
        fprintf('Generating model for the static state (Index %d)...\n', state_loop_idx);
        % current_state_vector remains empty, which is handled by createStateComponentList
    else
        error('InternalError:InvalidLoopStateForNoSwitches', 'Unexpected state_loop_idx %d when no switches are present and not processing static state 1.', state_loop_idx);
    end

    component_lists{state_loop_idx} = createStateComponentList(components, switches, current_state_vector, config, node_map_final, state_loop_idx);

    try
        [An, Bn, Xn_i, state_dim_n, mna_results_i] = generate_AB_for_state(component_lists{state_loop_idx}, node_map_final, coupling_map, config, u_vars, current_state_vector, state_loop_idx); 
    catch ME_AB
         fprintf(2, 'Error generating A/B matrices for state index %d (vector [%s]): %s\n', state_loop_idx, num2str(current_state_vector), ME_AB.message);
         fprintf(2, 'Occurred in function %s, line %d.\n', ME_AB.stack(1).name, ME_AB.stack(1).line);
         rethrow(ME_AB); 
    end

    % Consistency check for state variables, now based on the first *processed* state
    if k_proc_loop == 1 
        state_dim = state_dim_n;
        first_state_vars = Xn_i;
        if state_dim == 0 && length(indices_to_process) > 1 
             warning('AnalyzeSwitchingConverter:NoStateVarsFirstProcessed', ...
                'First processed state (index %d) resulted in zero state variables. Subsequent states might be inconsistent.', state_loop_idx);
        elseif state_dim == 0 && length(indices_to_process) == 1 && num_switches > 0 % Has switches, but chosen state has no dynamics
             fprintf('  The processed state (index %d) has no state variables (e.g., purely resistive for this configuration).\n', state_loop_idx);
        elseif state_dim == 0 && num_switches == 0 % No switches, no states
             fprintf('  Circuit has no state variables (static circuit).\n');
        end
    elseif state_dim_n ~= state_dim || ~isequal(sort(Xn_i), sort(first_state_vars))
        warning('AnalyzeSwitchingConverter:InconsistentStates', ...
                'State variables or dimension differ between processed switching states. State Index %d (%d vars): {%s} vs First Processed State (index %d, %d vars): {%s}', ...
                state_loop_idx, state_dim_n, strjoin(Xn_i,','), indices_to_process(1), state_dim, strjoin(first_state_vars,','));
        consistent_state_vars = false;
         error('State vector definition is inconsistent between the processed states. Cannot proceed.');
    end

    try
        [Cn, Dn] = derive_CD_for_state(config.outputs, component_lists{state_loop_idx}, mna_results_i, Xn_i, u_vars, config); 
    catch ME_CD
         fprintf(2, 'Error deriving C/D matrices for state index %d (vector [%s]): %s\n', state_loop_idx, num2str(current_state_vector), ME_CD.message);
         fprintf(2, 'Occurred in function %s, line %d.\n', ME_CD.stack(1).name, ME_CD.stack(1).line);
         rethrow(ME_CD); 
    end
    
    A_matrices{state_loop_idx} = An;
    B_matrices{state_loop_idx} = Bn;
    C_matrices{state_loop_idx} = Cn;
    D_matrices{state_loop_idx} = Dn;
    state_vars_list{state_loop_idx} = Xn_i; 

    fprintf('\n      --- State-Space Matrices for State Index %d ---\n', state_loop_idx);
    fprintf('      State Variables (Xn): %s\n', strjoin(Xn_i, ', '));
    fprintf('      A matrix (%d x %d):\n', size(An,1), size(An,2)); disp(An);
    fprintf('      B matrix (%d x %d):\n', size(Bn,1), size(Bn,2)); disp(Bn);
    fprintf('      C matrix (%d x %d):\n', size(Cn,1), size(Cn,2)); disp(Cn);
    fprintf('      D matrix (%d x %d):\n', size(Dn,1), size(Dn,2)); disp(Dn);
    fprintf('      --- End Matrices for State Index %d ---\n', state_loop_idx);

end % End loop through states to process

if ~consistent_state_vars && length(indices_to_process) > 1
    error('State vector definition is inconsistent between the processed states. Cannot proceed.');
end

% If any states were processed, use the state variables from the first one processed.
% If no states were processed (e.g., indices_to_process is empty), final_state_vars remains empty.
if ~isempty(indices_to_process) && ~isempty(first_state_vars)
    final_state_vars = first_state_vars; 
else
    final_state_vars = {}; % Ensure it's an empty cell if no states were processed or no state vars found
end


if ~isempty(indices_to_process)
    if state_dim > 0
        fprintf('State-space models generated for %d specified/valid states. State dimension: %d.\n', length(indices_to_process), state_dim);
        if ~isempty(final_state_vars)
            fprintf('State Variables: %s\n', strjoin(final_state_vars, ', '));
        else
             fprintf('State Variables: None identified from processed states.\n');
        end
    else % state_dim == 0
         fprintf('State-space models generated for %d specified/valid state(s), but circuit has no dynamic state variables.\n', length(indices_to_process));
    end
else
     fprintf('No states were processed for model generation (indices_to_process is empty).\n');
     final_state_vars = {}; % Ensure it's empty
     state_dim = 0; % No states processed means effective state dimension for further steps is 0
end




% =========================================================================
% --- Pre-calculate Models (Phase 1 & 2 done) ---
% The A, B, C, D matrices and final_state_vars are now available.
% =========================================================================

   % Initialize variables to be calculated inside the loop
    run_times_ms = [];
    timing_stats = struct();
    x0_ss = zeros(state_dim, 1); % Default initial condition
    
    % optimal_dt will store the per-switch deadtimes after optimization or fixed values.
    % It's initialized based on the config, which now holds per-switch values.
    optimal_dt = config.initial_deadtimes(:)'; % Ensure it's a row vector from per-switch config.initial_deadtimes
                                            % config.initial_deadtimes is set from user_deadtimes_applied_to_switches in runner
    
    dt_history = [];
    voltage_history = [];

    % Determine if running in timing mode
    is_timing_mode = config.num_timing_runs > 1;
    if is_timing_mode
        fprintf('\n--- Entering Timing Mode (%d runs) ---\n', config.num_timing_runs);
        run_times_ms = zeros(config.num_timing_runs, 1);
        if config.suppress_console_output
             fprintf('   (Console output during runs will be suppressed)\n');
        end
    end

% --- Main Calculation Loop (Runs once in normal mode, N times in timing mode) ---
    main_calc_timer = tic; 

    for run_idx = 1:config.num_timing_runs
        if is_timing_mode
            %fprintf('   Timing Run %d/%d...\n', run_idx, config.num_timing_runs);
            iter_timer = tic; 
        end

        % --- Create a local copy of params for this iteration ---
        params_iter = config; 
        
        % Add essential dynamic info to params_iter
        params_iter.state_dim = state_dim;
        params_iter.state_vars = final_state_vars; 
        params_iter.num_valid_states = num_valid_states;
        params_iter.num_switches = num_switches; 
        params_iter.u_vars = u_vars;
        params_iter.num_inputs = num_inputs;
        params_iter.parsed_input_values = parsed_input_vals;
        params_iter.valid_states = valid_states;
        
        % --- Phase 3: ESVA Setup & Deadtime Initialization (Inside Loop) ---
        if ~is_timing_mode || ~params_iter.suppress_console_output 
            fprintf('\n--- Phase 3: ESVA Setup & Deadtime Initialization ---\n');
        end


        topology_info = createTopologyInfo(valid_states, switches, final_state_vars, params_iter); 
        params_iter.topology_info = topology_info;
        
        params_iter = initializeDeadtimes(params_iter); 
        
        current_run_final_dts_to_use = params_iter.dt; 

        dt_history_iter = []; 
        voltage_history_iter = []; 

        % --- Phase 4: ZVS Deadtime Optimization (Conditional, Inside Loop) ---
        if params_iter.run_optimization
            if ~is_timing_mode || ~params_iter.suppress_console_output 
                fprintf('\n--- Phase 4: ZVS Deadtime Optimization ---\n');
                fprintf('Running ZVS Deadtime Optimization for per-switch deadtimes...\n');
            end
            [optimized_dts_from_nr, dt_hist_nr, voltage_hist_nr] = findOptimalDeadtimes(params_iter, A_matrices, B_matrices);
            
            current_run_final_dts_to_use = optimized_dts_from_nr; 
            dt_history_iter = dt_hist_nr; 
            voltage_history_iter = voltage_hist_nr;
            
            if ~is_timing_mode || ~params_iter.suppress_console_output
                fprintf('Optimization attempt complete for this run.\n');
                if ~isempty(current_run_final_dts_to_use) && ~any(isnan(current_run_final_dts_to_use))
                    fprintf('Per-switch deadtimes selected by optimizer: [%s] ns\n', sprintf('%.15e ', current_run_final_dts_to_use*1e9));
                else
                    fprintf('Optimizer did not return valid deadtimes. Check warnings from findOptimalDeadtimes.\n');
                end
            end
        else 
            if ~is_timing_mode || ~params_iter.suppress_console_output
                 fprintf('\n--- Phase 4: ZVS Deadtime Optimization ---\n');
                 fprintf('Skipping ZVS deadtime optimization. Using fixed per-switch deadtimes: [%s] ns\n', sprintf('%.15e ', current_run_final_dts_to_use*1e9));
            end
        end

        % % CRITICAL DEBUG: Log params_iter.dt before PSS calculation
        params_iter.dt = current_run_final_dts_to_use;
        params_iter.deadtimes_applied = current_run_final_dts_to_use; 
        % fprintf('DEBUG LOOP (run %d): params_iter.dt for PSS = [%s]\n', run_idx, sprintf('%.17e ', params_iter.dt));
        % assignin('base', sprintf('dbg_params_iter_dt_pss_run%d', run_idx), params_iter.dt);


        % --- Phase 5: Core PSS Calculation (Inside Loop) ---
        if ~is_timing_mode || ~params_iter.suppress_console_output 
             fprintf('\n--- Phase 5: Core PSS Calculation ---\n');
        end
        x0_ss_iter = zeros(state_dim, 1); 
        if state_dim > 0
            if ~is_timing_mode || ~params_iter.suppress_console_output 
                fprintf('Calculating steady-state initial conditions using ESVA...\n');
                fprintf('  Using applied per-switch deadtimes for sequence generation: [%s] ns\n', sprintf('%.17e ', params_iter.dt*1e9)); 
            end
            x0_ss_iter = enhancedSVA(params_iter, A_matrices, B_matrices); 

        % +++ START DIAGNOSTIC CALL TO simulateZVSConditions +++
        if ~params_iter.run_optimization && run_idx == 1 % Only for the non-optimization run, only once
            fprintf('\n+++ DIAGNOSTIC: Manually calling simulateZVSConditions with forced deadtimes +++\n');
            % Setup params specifically for this diagnostic call
            params_diag_zvs = params_iter; % Copy current params
            
            % Ensure is_switch_targeted_for_zvs and zvs_cap_state_var_indices are set up
            % This logic is similar to what's at the beginning of findOptimalDeadtimes
            num_sw_diag = params_diag_zvs.num_switches;
            is_targeted_diag = false(1, num_sw_diag);
            cap_indices_diag = NaN(1, num_sw_diag);

            if isfield(params_diag_zvs, 'topology_info') && ...
               isfield(params_diag_zvs.topology_info, 'zvs_target_map_switch_to_cap_idx') && ...
               isa(params_diag_zvs.topology_info.zvs_target_map_switch_to_cap_idx, 'containers.Map')
                
                map_sw_to_cap_diag = params_diag_zvs.topology_info.zvs_target_map_switch_to_cap_idx;
                
                % Determine which switches are targeted based on prms.zvs_switches_to_target
                % For this test, we assume M1 (switch 1) is the primary interest.
                % You might want to make this more general if needed, or ensure 
                % prms.zvs_switches_to_target is set in the runner even if optimization is off.
                % For simplicity, assume M1 is switch 1 and is targeted for this diagnostic.
                if num_sw_diag >= 1
                    % Manually set M1 as targeted for this diagnostic
                    % Find M1's index in switches_ordered_list
                    m1_idx_diag = 0;
                    for k_sw_find = 1:num_sw_diag
                        if strcmp(params_diag_zvs.switches_ordered_list(k_sw_find).name, 'M1')
                            m1_idx_diag = k_sw_find;
                            break;
                        end
                    end
                    
                    if m1_idx_diag > 0 % If M1 was found
                        is_targeted_diag(m1_idx_diag) = true; 
                        if isKey(map_sw_to_cap_diag, m1_idx_diag)
                            cap_indices_diag(m1_idx_diag) = map_sw_to_cap_diag(m1_idx_diag);
                        else
                            fprintf('DIAGNOSTIC: M1 not found in zvs_target_map_switch_to_cap_idx.\n');
                        end
                    else
                        fprintf('DIAGNOSTIC: Switch M1 not found in switches_ordered_list.\n');
                    end
                end
            else
                 fprintf('DIAGNOSTIC: Missing topology_info or zvs_target_map_switch_to_cap_idx for diagnostic.\n');
            end
            
            params_diag_zvs.is_switch_targeted_for_zvs = is_targeted_diag;
            params_diag_zvs.zvs_cap_state_var_indices = cap_indices_diag;
            params_diag_zvs.deadtimes_applied = params_iter.dt; % Use the forced deadtimes

            if any(is_targeted_diag) && ~any(isnan(cap_indices_diag(is_targeted_diag)))
                fprintf('DIAGNOSTIC: Calling simulateZVSConditions with deadtimes: [%s] ns\n', sprintf('%.3f ', params_diag_zvs.deadtimes_applied*1e9));
                fprintf('DIAGNOSTIC: x0_ss_iter for this call: [%s]\n', sprintf('%.6e ', x0_ss_iter));
                
                % Call simulateZVSConditions (ensure it has the fprintf statements for j_sw == 1)
                [diag_cap_voltages, diag_cap_rates] = simulateZVSConditions(params_diag_zvs, x0_ss_iter, A_matrices, B_matrices);
                
                fprintf('DIAGNOSTIC RESULTS from simulateZVSConditions:\n');
                if m1_idx_diag > 0 && m1_idx_diag <= length(diag_cap_voltages)
                    fprintf('  M1 (Switch %d) Forced ZVS Voltage: %.6e V\n', m1_idx_diag, diag_cap_voltages(m1_idx_diag));
                    fprintf('  M1 (Switch %d) Forced ZVS dV/dt:  %.6e V/s\n', m1_idx_diag, diag_cap_rates(m1_idx_diag));
                else
                     fprintf('  Could not retrieve M1 diagnostic results (index issue).\n');
                end
                assignin('base', 'dbg_diag_cap_voltages', diag_cap_voltages);
                assignin('base', 'dbg_diag_cap_rates', diag_cap_rates);
            else
                fprintf('DIAGNOSTIC: M1 not properly targeted or cap index missing for diagnostic call to simulateZVSConditions.\n');
            end
            fprintf('+++ END DIAGNOSTIC CALL +++\n\n');
        end
        % +++ END DIAGNOSTIC CALL TO simulateZVSConditions +++



            % % CRITICAL DEBUG: Log x0_ss_iter
            % assignin('base', sprintf('dbg_x0_ss_iter_run%d', run_idx), x0_ss_iter);
             % fprintf('DEBUG LOOP (run %d): x0_ss_iter = [%s]\n', run_idx, sprintf('%.17e ', x0_ss_iter));

            if ~is_timing_mode || ~params_iter.suppress_console_output 
                fprintf('Steady-state initial conditions calculated for this run.\n');
            end
        else
             if ~is_timing_mode || ~params_iter.suppress_console_output
                fprintf('Skipping PSS calculation as there are no state variables.\n');
             end
        end

        if run_idx == config.num_timing_runs
             x0_ss = x0_ss_iter;
             optimal_dt = current_run_final_dts_to_use; 
             
             if params_iter.run_optimization 
                dt_history = dt_history_iter;
                voltage_history = voltage_history_iter;
             else
                dt_history = []; 
                voltage_history = [];
             end
             params_final = params_iter; 
        end
        
        if is_timing_mode
            run_times_ms(run_idx) = toc(iter_timer) * 1000;
            if ~config.suppress_console_output 
                fprintf('   Run %d finished in %.4f ms\n', run_idx, run_times_ms(run_idx));
            end
        end

    end % End main calculation loop (for run_idx)

    total_calc_time_ms = toc(main_calc_timer) * 1000; 
    
    if is_timing_mode
        if ~config.suppress_console_output
            fprintf('--- Timing Mode Finished (Total calculation loop time: %.3f ms) ---\n', total_calc_time_ms);
        end
        timing_stats.num_runs = config.num_timing_runs;
        timing_stats.run_times_ms = run_times_ms;
        timing_stats.mean_time_ms = mean(run_times_ms);
        timing_stats.std_dev_time_ms = std(run_times_ms);
        timing_stats.min_time_ms = min(run_times_ms);
        timing_stats.max_time_ms = max(run_times_ms);
    else
        timing_stats = struct(); 
        if ~exist('params_final', 'var') 
            params_final = params_iter; 
        end
    end
    
    % CRITICAL DEBUG: Log optimal_dt and params_final.dt before waveform generation
    fprintf('DEBUG POST-LOOP: optimal_dt to be assigned to params_final.dt = [%s]\n', sprintf('%.17e ', optimal_dt));
    assignin('base', 'dbg_optimal_dt_post_loop', optimal_dt);
    
    params_final.dt = optimal_dt; 
    params_final.deadtimes_applied = optimal_dt;

    fprintf('DEBUG PRE-WAVEFORM: params_final.dt = [%s]\n', sprintf('%.17e ', params_final.dt));
    fprintf('DEBUG PRE-WAVEFORM: x0_ss for waveforms = [%s]\n', sprintf('%.17e ', x0_ss));
    assignin('base', 'dbg_params_final_dt_pre_waveform', params_final.dt);
    assignin('base', 'dbg_x0_ss_pre_waveform', x0_ss);
    
    % --- Phase 5: Output Generation (Run Once, after all calculation runs) ---
    if ~config.suppress_console_output 
        fprintf('\n--- Phase 5: Output Generation ---\n');
    end

    if ~config.suppress_console_output; fprintf('Generating steady-state waveforms (3 periods)...\n'); end
    
    [t_out, x_out, y_out] = generateOutputWaveforms(params_final, x0_ss, A_matrices, B_matrices, C_matrices, D_matrices); 
    if ~config.suppress_console_output; fprintf('Waveform generation complete.\n'); end

    time_before_packaging_ms = toc(timerVal) * 1000;

    results = packageResults(params_final, x0_ss, optimal_dt, dt_history, voltage_history, time_before_packaging_ms, t_out, x_out, y_out, valid_states, A_matrices, B_matrices, C_matrices, D_matrices, components, switches, node_map_final, timing_stats); 

    if ~config.suppress_console_output; displayResults(results); end
    if ~config.suppress_console_output; plotResults(results); end

    simulation_time_ms_total_func = toc(timerVal) * 1000; 
    results.simulation_time_ms = simulation_time_ms_total_func; 

    if ~config.suppress_console_output
        fprintf('\nAnalysis complete. Total function time (including plotting): %.3f ms\n', simulation_time_ms_total_func); 
        fprintf('=====================================================\n');
    end
    

    catch ME
        fprintf(2, '\n==================== ANALYSIS FAILED ====================\n');
        fprintf(2, 'Error Message: %s\n', ME.message);
        fprintf(2, '---------------------------------------------------------\n');
        fprintf(2, 'Error occurred in:\n');
        for k_err = 1:length(ME.stack)
            fprintf(2, '  File: %s\n', ME.stack(k_err).file);
            fprintf(2, '  Function: %s\n', ME.stack(k_err).name);
            fprintf(2, '  Line: %d\n', ME.stack(k_err).line);
            if k_err == 1 % Print loop index if error is in the main function's loop
                if strcmp(ME.stack(k_err).name, 'AnalyzeSwitchingConverter') && exist('state_loop_idx', 'var')
                     fprintf(2, '  (Occurred during processing state_loop_idx = %d)\n', state_loop_idx);
                end
            end
            fprintf(2, '---------------------------------------------------------\n');
        end
        fprintf(2, '=========================================================\n');
        % Return an empty or minimal results structure on failure
        results = struct('error', ME.message, 'error_stack', ME.stack);
    end % End main try-catch

% =========================================================================
% =========================================================================
%                            NESTED HELPER FUNCTIONS
% =========================================================================
% =========================================================================

 %% --- Phase 1 Helpers ---
function config_out = validateAndMergeConfig(config_in)
        % Validates required fields and merges with defaults.
        % MODIFIED: Added 'user_complementary_switch_pairs' to defaults and validation.
        
        fprintf('Validating and merging configuration...\n');
        config_out = config_in; % Start with user config

        % --- Define Default Configuration ---
        default_config = struct( ...
            'stop_after_phase1', false, ...
            'input_values', containers.Map(), ... 
            'include_parasitics', false, ...      
            'switch_parasitics', struct(), ...    
            'ts_switches', {{}}, ...              
            'deadtimes_applied', [], ...         
            'switches_ordered_list', struct('name', {}), ... 
            'user_operational_state_vectors', {{}}, ...
            'user_defined_state_index_sequence', [], ...
            'user_complementary_switch_pairs', {{}}, ...
            'time_tolerance_uniquetol', 1e-12, ... 
            'min_interval_duration', 1e-15, ...    
            'state_sequence', [], ... 
            'state_duration_proportions', containers.Map(), ... 
            'initial_deadtimes', [], ...        
            'run_optimization', false, ...
            'zvs_switches_to_target', {{}}, ...  
            'outputs', {{'all_states'}}, ...
            'max_iterations', 50, ... 
            'voltage_tolerance', 0.1, ...
            'damping_factor', 1.0, ...
            'min_dt', 1e-12, ...                   
            'max_dt_fraction', 0.2, ...            
            'max_individual_dt_fraction_of_T', 0.25, ... 
            'max_dt_step_change_fraction', 0.1, ... 
            'abstol', 1e-10, ...                   
            'degeneracy_threshold', 1e-20, ...     
            'significant_rate', 0, ...          
            'adjust_increase', 1.001, ...           
            'adjust_decrease', 0.999, ...           
            'num_timing_runs', 1, ... 
            'suppress_console_output', false);
      
        % --- Check Essential Physical Parameters (must be provided by runner) ---
        required_physical_fields = {'Ron', 'Roff', 'fs'};
        missing_req = {};
        for k_req = 1:length(required_physical_fields)
            field_name = required_physical_fields{k_req};
            if ~isfield(config_out, field_name) || isempty(config_out.(field_name)) || (isnumeric(config_out.(field_name)) && isnan(config_out.(field_name)))
                missing_req{end+1} = field_name;
            end
        end
        if ~isempty(missing_req)
             error('AnalyzeSwitchingConverter:MissingRequiredConfig', ...
                   'Config structure from runner script must include valid, non-empty fields: %s.', strjoin(missing_req, ', '));
        end
        if ~isscalar(config_out.Ron) || ~isnumeric(config_out.Ron) || config_out.Ron <= 0
            error('AnalyzeSwitchingConverter:InvalidRon', 'config.Ron must be a positive scalar number.');
        end
         if ~isscalar(config_out.Roff) || ~isnumeric(config_out.Roff) || config_out.Roff <= 0
            error('AnalyzeSwitchingConverter:InvalidRoff', 'config.Roff must be a positive scalar number.');
        end
         if ~isscalar(config_out.fs) || ~isnumeric(config_out.fs) || config_out.fs <= 0
            error('AnalyzeSwitchingConverter:InvalidFs', 'config.fs must be a positive scalar number.');
        end
        
        % Ensure T_period is consistent with fs, or set it if missing
        if ~isfield(config_out, 'T_period') || isempty(config_out.T_period) 
            config_out.T_period = 1/config_out.fs;
        elseif abs(config_out.T_period - 1/config_out.fs) > (1/config_out.fs * 1e-9) % Check with tolerance
             warning('AnalyzeSwitchingConverter:TperiodFsMismatch', 'config.T_period (%.4e) from runner script does not match 1/config.fs (%.4e). Using 1/fs.', config_out.T_period, 1/config_out.fs);
             config_out.T_period = 1/config_out.fs;
        end

        % --- Merge User Config with Defaults ---
        config_fields = fieldnames(default_config);
        allow_empty_fields = {'state_sequence', 'initial_deadtimes', ...
                              'outputs', 'zvs_switches_to_target', ...
                              'switch_parasitics', 'ts_switches', ...
                              'deadtimes_applied', 'switches_ordered_list', ...
                              'input_values', 'state_duration_proportions', ...
                              'user_operational_state_vectors', ...
                              'user_defined_state_index_sequence', ...
                              'user_complementary_switch_pairs'};
                          
        for i_cfg = 1:length(config_fields)
            field = config_fields{i_cfg};
            if ~isfield(config_out, field)
                config_out.(field) = default_config.(field); 
            elseif isempty(config_out.(field)) && ~isstruct(config_out.(field)) && ~isa(config_out.(field), 'containers.Map')
                % Allow specific fields to be empty, otherwise use default
                if ~ismember(field, allow_empty_fields)
                    config_out.(field) = default_config.(field); 
                end
            end
        end
        
        % --- Ensure specific parameters from common_params in runner script override defaults if they exist in config_in ---
        if isfield(config_in, 'time_tolerance_uniquetol') && ~isempty(config_in.time_tolerance_uniquetol)
            config_out.time_tolerance_uniquetol = config_in.time_tolerance_uniquetol;
        end
        if isfield(config_in, 'min_interval_duration') && ~isempty(config_in.min_interval_duration)
            config_out.min_interval_duration = config_in.min_interval_duration;
        end
        if isfield(config_in, 'max_individual_dt_fraction_of_T') && ~isempty(config_in.max_individual_dt_fraction_of_T)
            config_out.max_individual_dt_fraction_of_T = config_in.max_individual_dt_fraction_of_T;
        end
         if isfield(config_in, 'max_dt_step_change_fraction') && ~isempty(config_in.max_dt_step_change_fraction)
            config_out.max_dt_step_change_fraction = config_in.max_dt_step_change_fraction;
        end

        % --- Validate Specific Config Fields (Basic Types and General Logic) ---
        if ~isa(config_out.input_values, 'containers.Map')
            error('AnalyzeSwitchingConverter:InvalidInputValuesType', 'config.input_values must be a containers.Map.');
        end
        
        if ~islogical(config_out.include_parasitics) || ~isscalar(config_out.include_parasitics)
             error('AnalyzeSwitchingConverter:InvalidIncludeParasitics', 'config.include_parasitics must be a logical scalar.');
        end
        if config_out.include_parasitics
             if ~isfield(config_out, 'switch_parasitics') || ~isstruct(config_out.switch_parasitics)
                 error('AnalyzeSwitchingConverter:MissingSwitchParasiticsStruct', 'If config.include_parasitics is true, config.switch_parasitics must be provided as a struct from the runner script.');
             end
        end

        if ~iscell(config_out.user_operational_state_vectors)
            error('AnalyzeSwitchingConverter:InvalidUserOpStatesType', 'config.user_operational_state_vectors must be a cell array.');
        end
        
        if ~isnumeric(config_out.user_defined_state_index_sequence) || (~isempty(config_out.user_defined_state_index_sequence) && ~isvector(config_out.user_defined_state_index_sequence))
            error('AnalyzeSwitchingConverter:InvalidUserDefSeqType', 'config.user_defined_state_index_sequence must be a numeric vector (or empty).');
        end
        if ~isempty(config_out.user_defined_state_index_sequence) && any(config_out.user_defined_state_index_sequence < 1 | mod(config_out.user_defined_state_index_sequence, 1) ~= 0)
            error('AnalyzeSwitchingConverter:InvalidUserDefSeqContent', 'Elements of config.user_defined_state_index_sequence must be positive integers.');
        end

        % ADDED VALIDATION for user_complementary_switch_pairs
        if ~iscell(config_out.user_complementary_switch_pairs)
             error('AnalyzeSwitchingConverter:InvalidUserCompPairsType', 'config.user_complementary_switch_pairs must be a cell array (e.g., {{''M1'',''M2''},{''M3'',''M4''}}).');
        end
        for k_pair_val = 1:length(config_out.user_complementary_switch_pairs)
            pair_val = config_out.user_complementary_switch_pairs{k_pair_val};
            if ~iscell(pair_val) || length(pair_val) ~= 2 || ~ischar(pair_val{1}) || ~ischar(pair_val{2})
                 error('AnalyzeSwitchingConverter:InvalidUserCompPairFormat', 'Each element in config.user_complementary_switch_pairs must be a cell array of two switch name strings.');
            end
        end
        % END ADDED VALIDATION

        if ~isnumeric(config_out.num_timing_runs) || ~isscalar(config_out.num_timing_runs) || config_out.num_timing_runs < 1 || mod(config_out.num_timing_runs, 1) ~= 0
            error('AnalyzeSwitchingConverter:InvalidNumTimingRuns', 'config.num_timing_runs must be a positive integer scalar.');
        end
        if ~islogical(config_out.suppress_console_output) || ~isscalar(config_out.suppress_console_output)
            error('AnalyzeSwitchingConverter:InvalidSuppressOutput', 'config.suppress_console_output must be a logical scalar (true or false).');
        end
        
        if ~isscalar(config_out.time_tolerance_uniquetol) || ~isnumeric(config_out.time_tolerance_uniquetol) || config_out.time_tolerance_uniquetol < 0 
            error('AnalyzeSwitchingConverter:InvalidTimeToleranceUniquetol', 'config.time_tolerance_uniquetol must be a non-negative scalar number.');
        end
        if ~isscalar(config_out.min_interval_duration) || ~isnumeric(config_out.min_interval_duration) || config_out.min_interval_duration < 0 
            error('AnalyzeSwitchingConverter:InvalidMinIntervalDuration', 'config.min_interval_duration must be a non-negative scalar number.');
        end

        % --- Type checks for switch-dependent parameters (Length checks deferred) ---
        if ~config_out.stop_after_phase1
            if ~config_out.suppress_console_output
                fprintf('  Performing full configuration validation (stop_after_phase1 is false).\n');
                fprintf('  (Note: Length checks for switch-dependent arrays like ts_switches, deadtimes_applied will occur after parseNetlist).\n');
            end

            % Debug print for ts_switches
            if ~config_out.suppress_console_output
                fprintf('  DEBUG validateAndMergeConfig: Before ts_switches type check. Class: %s, IsEmpty: %d\n', class(config_out.ts_switches), isempty(config_out.ts_switches));
            end
            if ~iscell(config_out.ts_switches) % This should be fine even if empty {}
                error('AnalyzeSwitchingConverter:InvalidTsSwitchesType', 'config.ts_switches must be a cell array.');
            end
            % Individual cell content validation
            for i_ts_val = 1:length(config_out.ts_switches) % Loop is fine if length is 0
                if ~isnumeric(config_out.ts_switches{i_ts_val}) || ~isvector(config_out.ts_switches{i_ts_val})
                     error('AnalyzeSwitchingConverter:InvalidTsSwitchesCellContent', ...
                          'Each cell of config.ts_switches (e.g., for switch %d) must be a numeric vector.', i_ts_val);
                end
                time_tol_check_val = config_out.T_period * 1e-12; 
                if any(config_out.ts_switches{i_ts_val} < -time_tol_check_val) || any(config_out.ts_switches{i_ts_val} > config_out.T_period + time_tol_check_val) 
                     error('AnalyzeSwitchingConverter:InvalidIdealTimeRangeInConfig', ...
                          'Ideal switching instants in config.ts_switches{%d} (value: [%s]) must be within [0, T_period=%.4es].', ...
                          i_ts_val, num2str(config_out.ts_switches{i_ts_val}), config_out.T_period);
                end
            end
            
            % Debug print for deadtimes_applied
            if ~config_out.suppress_console_output
                 fprintf('  DEBUG validateAndMergeConfig: Before deadtimes_applied type check. Class: %s, IsEmpty: %d, Value: %s\n', class(config_out.deadtimes_applied), isempty(config_out.deadtimes_applied), mat2str(config_out.deadtimes_applied));
            end
            % Corrected check: Allow empty numeric vector
            if ~isempty(config_out.deadtimes_applied) && (~isnumeric(config_out.deadtimes_applied) || ~isvector(config_out.deadtimes_applied))
                error('AnalyzeSwitchingConverter:InvalidDeadtimesAppliedType', 'If config.deadtimes_applied is not empty, it must be a numeric vector.');
            elseif isempty(config_out.deadtimes_applied) && ~isnumeric(config_out.deadtimes_applied) % Handles cases like non-numeric empty
                 error('AnalyzeSwitchingConverter:InvalidEmptyDeadtimesAppliedType', 'config.deadtimes_applied must be a numeric vector (or an empty numeric array []).');
            end
            if ~isempty(config_out.deadtimes_applied) && any(config_out.deadtimes_applied < 0) 
                error('AnalyzeSwitchingConverter:NegativeDeadtimeApplied', 'Values in config.deadtimes_applied must be non-negative.');
            end
            
            if ~isstruct(config_out.switches_ordered_list) || (~isempty(config_out.switches_ordered_list) && ~isfield(config_out.switches_ordered_list, 'name'))
                 error('AnalyzeSwitchingConverter:InvalidSwitchesOrderedList', 'config.switches_ordered_list must be a struct array with a "name" field, or an empty struct array.');
            end
            
            % Corrected check: Allow empty numeric vector for initial_deadtimes
            if ~isempty(config_out.initial_deadtimes) && (~isnumeric(config_out.initial_deadtimes) || ~isvector(config_out.initial_deadtimes))
                error('AnalyzeSwitchingConverter:InvalidInitialDTFormat', 'If config.initial_deadtimes is not empty, it must be a numeric vector.');
            elseif isempty(config_out.initial_deadtimes) && ~isnumeric(config_out.initial_deadtimes)
                 error('AnalyzeSwitchingConverter:InvalidEmptyInitialDTFormat', 'config.initial_deadtimes must be a numeric vector (or an empty numeric array []).');
            end
            if ~isempty(config_out.initial_deadtimes) && any(config_out.initial_deadtimes < 0)
                 error('AnalyzeSwitchingConverter:NegativeInitialDeadtime', 'Values in config.initial_deadtimes must be non-negative.');
            end
            
            % ... (other type validations that don't depend on num_switches)
            scalar_positive_fields = {'max_iterations', 'voltage_tolerance', 'damping_factor', ...
                                      'max_dt_fraction', 'abstol', 'degeneracy_threshold', ...
                                      'adjust_increase', 'adjust_decrease', ...
                                      'max_individual_dt_fraction_of_T', 'max_dt_step_change_fraction'};
            scalar_nonnegative_fields = {'min_dt', 'significant_rate'};
            
            for f_idx = 1:length(scalar_positive_fields)
                fname = scalar_positive_fields{f_idx}; 
                if ~isscalar(config_out.(fname)) || ~isnumeric(config_out.(fname)) || config_out.(fname) <= 0
                    error('AnalyzeSwitchingConverter:InvalidNumericParamPositive', 'config.%s must be a positive scalar number. Value: %g', fname, config_out.(fname)); 
                end
            end
            for f_idx = 1:length(scalar_nonnegative_fields)
                fname = scalar_nonnegative_fields{f_idx}; 
                if ~isscalar(config_out.(fname)) || ~isnumeric(config_out.(fname)) || config_out.(fname) < 0
                    error('AnalyzeSwitchingConverter:InvalidNumericParamNonNegative', 'config.%s must be a non-negative scalar number. Value: %g', fname, config_out.(fname)); 
                end
            end
            if config_out.max_individual_dt_fraction_of_T >= 1.0
                warning('AnalyzeSwitchingConverter:LargeMaxIndividualDT', 'config.max_individual_dt_fraction_of_T (%.2f) is >= 1.0. This means a single deadtime can be >= T_period, which is unusual.', config_out.max_individual_dt_fraction_of_T);
            end
             if config_out.max_dt_step_change_fraction >= 1.0
                warning('AnalyzeSwitchingConverter:LargeMaxStepChangeDT', 'config.max_dt_step_change_fraction (%.2f) is >= 1.0. This means a single deadtime can change by >= T_period in one ZVS iteration.', config_out.max_dt_step_change_fraction);
            end

        else 
             if ~config_out.suppress_console_output
                fprintf('  Skipping full configuration validation (stop_after_phase1 is true).\n');
             end
        end 
        
        if ~config_out.suppress_console_output
            fprintf('Configuration validated and merged.\n');
        end
    end % end validateAndMergeConfig

 % --- Validate lengths of per-switch config fields against num_switches ---
        if ~config.stop_after_phase1 % Only for full run
            
            % This block now relies on 'num_switches' being correctly defined and in scope
            % from the definition block that follows the 'parseNetlist' call.
            
            fprintf('DEBUG VALIDATION BLOCK (Simplified): config.stop_after_phase1 is false.\n');
            % Check if num_switches itself exists, as a basic sanity check for this block
            if ~exist('num_switches', 'var')
                error('AnalyzeSwitchingConverter:NumSwitchesMissingForValidation', "'num_switches' variable not found. It should have been defined after parseNetlist.");
            end
            fprintf('DEBUG VALIDATION BLOCK (Simplified): Using num_switches = %d for length checks.\n', num_switches);
            
            if ~config.suppress_console_output && num_switches > 0 
                fprintf('  Validating per-switch configuration field lengths against %d identified switches.\n', num_switches);
            elseif ~config.suppress_console_output && num_switches == 0
                 fprintf('  Skipping per-switch config length validation as num_switches is 0.\n');
            end
    
            % Validate config.ts_switches length
            if length(config.ts_switches) ~= num_switches
                error('AnalyzeSwitchingConverter:TsSwitchesLengthMismatch', ...
                      'Length of config.ts_switches (%d) must match the number of identified switches (%d).', ...
                      length(config.ts_switches), num_switches);
            end
    
            % Validate config.deadtimes_applied length
            if length(config.deadtimes_applied) ~= num_switches
                error('AnalyzeSwitchingConverter:DeadtimesAppliedLengthMismatch', ...
                      'Length of config.deadtimes_applied (%d) must match the number of identified switches (%d).', ...
                      length(config.deadtimes_applied), num_switches);
            end
            
            % Validate config.initial_deadtimes length (as it's now per-switch)
            if length(config.initial_deadtimes) ~= num_switches
                error('AnalyzeSwitchingConverter:InitialDeadtimesLengthMismatch', ...
                      'Length of config.initial_deadtimes (%d) must match the number of identified switches (%d).', ...
                      length(config.initial_deadtimes), num_switches);
            end
    
            % Validate config.switches_ordered_list length
            if length(config.switches_ordered_list) ~= num_switches
                error('AnalyzeSwitchingConverter:SwitchesOrderedListLengthMismatch', ...
                      'Length of config.switches_ordered_list (%d) must match the number of identified switches (%d).', ...
                      length(config.switches_ordered_list), num_switches);
            end
        end % end if ~config.stop_after_phase1
    %------------------------------------------------
    function [components_out, switches_out, nodeMap_out, coupling_map_out, config_updated, u_vars_out, parsed_input_values] = parseNetlist(filePath, config_in)
        % Parses netlist file (R, L, C, V, I, E, G, H, F, K, M).
        % Identifies components, switches, nodes, coupling, and input sources.
        % Stores parsed values for V/I sources.
        % Errors out if no independent V or I source is found.

        % --- Initialization ---
        components_out = {};
        switches_out = struct('name', {}, 'type', {}, 'nodes', {}, 'ctrl_node', {}, 'index', {});
        node_list = {'0'}; % Ground node '0' is always index 0
        nodeMap_out = containers.Map('KeyType', 'char', 'ValueType', 'double');
        nodeMap_out('0') = 0;
        node_idx_counter = 0; % Counter for non-ground nodes (1-based index)
        comp_idx = 0; % Counter for components parsed
        coupling_map_out = containers.Map('KeyType', 'char', 'ValueType', 'any'); % Stores {L1_name, L2_name, k} for coupled inductors
        config_updated = config_in; % Pass config through (no modifications here currently)
        u_vars_out = {}; % Cell array to store names of independent V and I sources found
        parsed_input_values = containers.Map('KeyType', 'char', 'ValueType', 'double'); % Stores parsed values for V/I sources

        % --- Read File ---
        try
            netlistCell = readlines(filePath);
            netlistCell = strtrim(netlistCell);
            % Remove empty lines and standard SPICE comment lines (*, ;)
            netlistCell = netlistCell(strlength(netlistCell) > 0);
            netlistCell = netlistCell(~startsWith(netlistCell, {'*', ';'}));
             % Remove other common comment styles (#, .) and MATLAB comments (%)
            netlistCell = netlistCell(~startsWith(netlistCell, '#'));
            netlistCell = netlistCell(~startsWith(netlistCell, '.'));
            netlistCell = netlistCell(~startsWith(netlistCell, '%')); % Added check for MATLAB comments

            if isempty(netlistCell)
                 error('AnalyzeSwitchingConverter:NetlistFileEmpty', 'Netlist file is empty or contains no valid component lines: %s', filePath);
            end
            fprintf('  Successfully read %d lines from %s.\n', length(netlistCell), filePath);
        catch ME
            error('AnalyzeSwitchingConverter:NetlistFileReadError', 'Failed to read netlist file "%s": %s', filePath, ME.message);
        end

        % --- Parse Lines ---
        % Iterate through each line of the processed netlist content.
        for line_idx = 1:length(netlistCell)
            line = strtrim(netlistCell{line_idx});
            parts = strsplit(line); % Split line into parts based on whitespace

            % Basic check for minimum required parts (e.g., K L1 L2 k needs 4)
            if length(parts) < 3
                warning('Skipping invalid netlist line %d (too few parts): %s', line_idx, line);
                continue;
            end

            comp_idx = comp_idx + 1; % Increment component counter
            comp = struct(); % Initialize component structure
            comp.name = parts{1};
            comp.type = upper(comp.name(1)); % First letter determines type

            % Validate component type character
            if ~isletter(comp.type)
                 warning('Skipping invalid netlist line %d ("%s") due to error: Invalid component type character "%s"', line_idx, line, comp.type);
                 comp_idx = comp_idx - 1; % Decrement counter as component is skipped
                 continue;
            end

            % Initialize flags and fields for the component struct
            comp.index = comp_idx; % Store original parse order index
            comp.is_resistor = false; comp.is_inductor = false; comp.is_capacitor = false;
            comp.is_voltage_source = false; comp.is_current_source = false;
            comp.is_vcvs = false; comp.is_vccs = false; comp.is_ccvs = false;
            comp.is_cccs = false; comp.is_mutual_ind = false; comp.is_switch = false;
            comp.nodes = {}; comp.controlling_nodes = {}; comp.controlling_element = '';
            comp.value = NaN; % Holds R, L, C value, V/I value, Gain/Gm/R/Beta, or k coeff
            current_nodes_to_map = {}; % Nodes associated with this component

            try % Use try-catch for parsing of each line
                % --- Parse based on component type ---
                switch comp.type
                    case {'R', 'L', 'C'} % Passive Components
                        if length(parts) < 4, error('Format Error: Needs Name N1 N2 Value'); end
                        comp.nodes = {parts{2}, parts{3}};
                        comp.value = parse_value_simple(parts{4}); % Use helper to parse value with suffixes
                        if isnan(comp.value), error('Value Error: Invalid numeric value'); end
                        % Set type flag and validate value
                        if comp.type == 'R', comp.is_resistor = true; if comp.value <= 0, error('Value Error: R must be > 0'); end; end
                        if comp.type == 'L', comp.is_inductor = true; if comp.value <= 0, error('Value Error: L must be > 0'); end; end
                        if comp.type == 'C', comp.is_capacitor = true; if comp.value <= 0, error('Value Error: C must be > 0'); end; end
                        current_nodes_to_map = comp.nodes;

                    case 'V' % Independent Voltage Source
                        if length(parts) < 4, error('Format Error: Needs Name N+ N- Value'); end
                        comp.nodes = {parts{2}, parts{3}};
                        comp.value = parse_value_simple(parts{4});
                        if isnan(comp.value), error('Value Error: Invalid numeric value'); end
                        comp.is_voltage_source = true;
                        current_nodes_to_map = comp.nodes;
                        u_vars_out{end+1} = comp.name; % Add name to list of inputs
                        parsed_input_values(comp.name) = comp.value; % Store parsed value

                    case 'I' % Independent Current Source
                        if length(parts) < 4, error('Format Error: Needs Name N+ N- Value'); end
                        comp.nodes = {parts{2}, parts{3}};
                        comp.value = parse_value_simple(parts{4});
                        if isnan(comp.value), error('Value Error: Invalid numeric value'); end
                        comp.is_current_source = true;
                        current_nodes_to_map = comp.nodes;
                        u_vars_out{end+1} = comp.name; % Add name to list of inputs
                        parsed_input_values(comp.name) = comp.value; % Store parsed value

                    case 'M' % Switch (MOSFET)
                        % Expects SPICE-like format: Mname Drain Gate Source Model
                        % Bulk node is ignored if present.
                        if length(parts) < 5, error('Format Error: MOSFET needs Name D G S Model'); end
                        drain_node = parts{2}; gate_node = parts{3}; source_node = parts{4};
                        comp.type = 'SW'; % Internal representation type
                        comp.is_switch = true;
                        comp.nodes = {drain_node, source_node}; % Store D-S nodes
                        comp.value = NaN; % Value not used directly for switch
                        comp.controlling_nodes = {gate_node}; % Store gate node
                        current_nodes_to_map = {drain_node, gate_node, source_node}; % All nodes need mapping
                        % Store switch-specific info
                        sw_idx = length(switches_out) + 1;
                        switches_out(sw_idx).name = comp.name;
                        switches_out(sw_idx).type = 'MOSFET'; % Could be extended later
                        switches_out(sw_idx).nodes = {drain_node, gate_node, source_node}; % D, G, S
                        switches_out(sw_idx).ctrl_node = gate_node; % Explicit control node
                        switches_out(sw_idx).index = comp_idx; % Index in the main component list

                    case 'E' % Voltage-Controlled Voltage Source (VCVS)
                        if length(parts) < 6, error('Format Error: Needs Name N+ N- NC+ NC- Gain'); end
                        comp.nodes = {parts{2}, parts{3}}; % Output nodes
                        comp.controlling_nodes = {parts{4}, parts{5}}; % Controlling voltage nodes
                        comp.value = parse_value_simple(parts{6}); % Gain value
                        if isnan(comp.value), error('Value Error: Invalid numeric gain'); end
                        comp.is_vcvs = true;
                        comp.is_voltage_source = true; % Treat as voltage source for MNA current variable
                        current_nodes_to_map = [comp.nodes, comp.controlling_nodes];

                    case 'G' % Voltage-Controlled Current Source (VCCS)
                        if length(parts) < 6, error('Format Error: Needs Name N+ N- NC+ NC- Gm'); end
                        comp.nodes = {parts{2}, parts{3}}; % Output nodes
                        comp.controlling_nodes = {parts{4}, parts{5}}; % Controlling voltage nodes
                        comp.value = parse_value_simple(parts{6}); % Transconductance Gm
                        if isnan(comp.value), error('Value Error: Invalid numeric Gm'); end
                        comp.is_vccs = true;
                        comp.is_current_source = true; % Treat as current source for MNA stamping
                        current_nodes_to_map = [comp.nodes, comp.controlling_nodes];

                    case 'H' % Current-Controlled Voltage Source (CCVS)
                        if length(parts) < 5, error('Format Error: Needs Name N+ N- Vctrl R'); end
                        comp.nodes = {parts{2}, parts{3}}; % Output nodes
                        comp.controlling_element = parts{4}; % Name of V, E, H source providing control current
                        comp.value = parse_value_simple(parts{5}); % Transresistance R
                        if isnan(comp.value), error('Value Error: Invalid numeric R'); end
                        comp.is_ccvs = true;
                        comp.is_voltage_source = true; % Treat as voltage source for MNA current variable
                        current_nodes_to_map = comp.nodes; % Only output nodes need mapping here

                    case 'F' % Current-Controlled Current Source (CCCS)
                        if length(parts) < 5, error('Format Error: Needs Name N+ N- Vctrl Beta'); end
                        comp.nodes = {parts{2}, parts{3}}; % Output nodes
                        comp.controlling_element = parts{4}; % Name of V, E, H source providing control current
                        comp.value = parse_value_simple(parts{5}); % Current gain Beta
                        if isnan(comp.value), error('Value Error: Invalid numeric Beta'); end
                        comp.is_cccs = true;
                        comp.is_current_source = true; % Treat as current source for MNA stamping
                        current_nodes_to_map = comp.nodes; % Only output nodes need mapping here

                    case 'K' % Mutual Inductance
                        if length(parts) < 4, error('Format Error: Needs Name L1 L2 k'); end
                        k_info = struct(); % Create struct to hold coupling info
                        k_info.L1_name = parts{2};
                        k_info.L2_name = parts{3};
                        k_info.k = parse_value_simple(parts{4}); % Coupling coefficient k
                        if isnan(k_info.k), error('Value Error: Invalid numeric k'); end
                        % Validate coupling coefficient range
                        if k_info.k < 0 || k_info.k > 1
                           warning('AnalyzeSwitchingConverter:InvalidCouplingCoeff', 'Coupling coefficient k=%.3f for %s must be [0, 1]. Clamping.', k_info.k, comp.name);
                           k_info.k = max(0, min(1, k_info.k)); % Clamp k to [0, 1]
                        end
                        comp.is_mutual_ind = true;
                        comp.value = k_info.k; % Store k value in component struct
                        comp.nodes = {k_info.L1_name, k_info.L2_name}; % Store coupled inductor names
                        % Add to coupling map for easy lookup later, mapping both ways
                        if isKey(coupling_map_out, k_info.L1_name) || isKey(coupling_map_out, k_info.L2_name)
                           warning('AnalyzeSwitchingConverter:DuplicateCoupling', 'Coupling involving %s or %s already defined. Overwriting previous definition with %s.', k_info.L1_name, k_info.L2_name, comp.name);
                        end
                        coupling_map_out(k_info.L1_name) = k_info; % Map L1 -> {L1, L2, k}
                        coupling_map_out(k_info.L2_name) = k_info; % Map L2 -> {L1, L2, k}
                        current_nodes_to_map = {}; % Nodes are inductor names, not mapped here

                    otherwise % Unknown component type
                        error('Unknown Type: Component name starts with "%s"', comp.type);
                end

                % --- Add nodes to map and list ---
                % Process nodes associated with the current component (excluding K elements)
                for j = 1:length(current_nodes_to_map)
                    node_name = current_nodes_to_map{j};

                    % Basic check for valid node names (alphanumeric/underscore or numeric string)
                    % Allows '0' for ground. Issues warning for non-standard names.
                    isNumericNode = ~isnan(str2double(node_name));
                    if ~strcmp(node_name, '0') && ~isvarname(node_name) && ~isNumericNode
                         warning('AnalyzeSwitchingConverter:NodeNameInvalid', 'Node name "%s" for component %s is not standard (alphanumeric/underscore or numeric). Check netlist.', node_name, comp.name);
                    end

                    % Add node to map if it's new and not ground
                    if ~isKey(nodeMap_out, node_name)
                        if ~strcmp(node_name, '0')
                            node_idx_counter = node_idx_counter + 1; % Increment 1-based index
                            nodeMap_out(node_name) = node_idx_counter; % Map name to index
                            node_list{end+1} = node_name; % Add to list of unique node names
                        end
                        % Ground node '0' is already mapped to index 0
                    end
                end
                components_out{end+1} = comp; % Add successfully parsed component to list

            catch ME % Catch parsing errors for the current line
                comp_idx = comp_idx - 1; % Decrement counter as component parsing failed
                warning('AnalyzeSwitchingConverter:ParseError', 'Skipping netlist line %d ("%s") due to error: %s', line_idx, line, ME.message);
                 % Attempt to remove from coupling map if K parsing failed
                 if comp.is_mutual_ind && exist('k_info','var')
                     if isfield(k_info,'L1_name') && isKey(coupling_map_out, k_info.L1_name); remove(coupling_map_out, k_info.L1_name); end
                      if isfield(k_info,'L2_name') && isKey(coupling_map_out, k_info.L2_name); remove(coupling_map_out, k_info.L2_name); end
                 end
                continue; % Skip to next line
            end
        end % End loop through netlist lines

        % --- Validate K element references ---
        % Ensure inductors mentioned in K lines exist in the component list
        k_keys = keys(coupling_map_out);
        inductor_names = {};
        for k_comp=1:length(components_out); if components_out{k_comp}.is_inductor; inductor_names{end+1} = components_out{k_comp}.name; end; end
        keys_to_remove = {}; processed_pairs = containers.Map; % Track processed pairs to avoid double warnings/removals
        for k_key=1:length(k_keys)
             key = k_keys{k_key};
             % Check if key still exists (might have been removed by previous iteration)
             if ~isKey(coupling_map_out, key), continue; end
             k_info = coupling_map_out(key);
             pair_key = strjoin(sort({k_info.L1_name, k_info.L2_name}), '-'); % Unique key for the pair

             % Process each pair only once
             if ~isKey(processed_pairs, pair_key)
                 processed_pairs(pair_key) = true; % Mark pair as processed
                 % Check if both inductors exist
                 if ~ismember(k_info.L1_name, inductor_names) || ~ismember(k_info.L2_name, inductor_names)
                     warning('AnalyzeSwitchingConverter:InvalidKElement', 'K element coupling %s and %s is invalid because one or both inductors were not found or parsed correctly. Removing coupling.', k_info.L1_name, k_info.L2_name);
                     % Mark both keys associated with this invalid pair for removal
                     keys_to_remove{end+1} = k_info.L1_name;
                     keys_to_remove{end+1} = k_info.L2_name;
                 end
             end
        end
        % Remove invalid coupling entries from the map
        keys_to_remove = unique(keys_to_remove);
        for k_rem = 1:length(keys_to_remove)
            if isKey(coupling_map_out, keys_to_remove{k_rem}); remove(coupling_map_out, keys_to_remove{k_rem}); end
        end

        % --- Finalize Input Source Check ---
        % Check if any independent V or I sources were found. Error if none.
        if isempty(u_vars_out)
             error('AnalyzeSwitchingConverter:NoInputSource', 'Parsing failed to find any independent voltage (V) or current (I) sources in the netlist. At least one independent source is required for analysis.');
        else
            fprintf('  Found %d independent source(s): %s\n', length(u_vars_out), strjoin(u_vars_out, ', '));
        end

    end % end parseNetlist

    %------------------------------------------------
 function [valid_s, num_valid] = generateValidStates(switches_in, cfg)
        % Generates all 2^N switch states.
        % Filters out the all-ON state ONLY IF n_sw >= 2, as an assumed
        % shoot-through condition. For n_sw = 1, both [0] and [1] are valid.

        n_sw = length(switches_in); % Number of switches identified
        if n_sw == 0
            % If there are no switches, the only "state" is the static circuit.
            % valid_s should be empty, and num_valid should represent the single static state.
            % The MNA for a static circuit will be handled by generate_AB_for_state.
            % The "state vector" for a no-switch case can be considered empty or a placeholder.
            % The actual "state index" for this single static state will be 1.
            valid_s = []; % No switch state vectors per se
            num_valid = 1;  % Represents the single static operating state
            if ~cfg.suppress_console_output
                fprintf('  No switches found. Circuit has 1 static operating state.\n');
            end
            return;
        end

        num_possible_states = 2^n_sw; % Total possible ON/OFF combinations
        % Generate all combinations (rows are states, columns are switches)
        % 'left-msb' ensures the first column corresponds to the first switch
        all_states = de2bi(0:(num_possible_states-1), n_sw, 'left-msb');

        valid_mask = true(num_possible_states, 1); % Start assuming all states are valid
        num_removed = 0; % Counter for removed states

        % --- Automatically remove the all-ON state ONLY IF n_sw >= 2 ---
        if n_sw >= 2 % Only consider shoot-through for 2 or more switches
            all_on_state_idx = find(all(all_states == 1, 2)); % Find row where all elements are 1

            if ~isempty(all_on_state_idx)
                % Check if this state hasn't already been invalidated (shouldn't happen here)
                if valid_mask(all_on_state_idx)
                    if ~cfg.suppress_console_output
                        fprintf('  Auto-removing all-ON state [%s] as shoot-through (n_sw >= 2).\n', num2str(ones(1,n_sw)));
                    end
                    valid_mask(all_on_state_idx) = false; % Mark the all-ON state as invalid
                    num_removed = num_removed + 1;
                end
            else
                % This case should not be reached if n_sw >= 2, as an all-ON state will always exist.
                if ~cfg.suppress_console_output
                     fprintf('  All-ON state not found among possible states (this is unusual for n_sw >= 2).\n');
                end
            end
        elseif ~cfg.suppress_console_output && n_sw == 1
            fprintf('  Single switch found. Both ON [1] and OFF [0] states are considered valid by default.\n');
        end
        % --- End Automatic Removal ---

        % --- Finalize valid states ---
        valid_s = all_states(valid_mask, :); % Keep only rows marked as valid
        
        % Ensure states are sorted if that's an implicit assumption elsewhere,
        % though ismember(..., 'rows') doesn't require it.
        % For consistency, sort by the decimal equivalent of the binary rows.
        if ~isempty(valid_s)
            [~, sort_order] = sortrows(valid_s); % Sort based on the rows themselves
            valid_s = valid_s(sort_order, :);
        end
        
        num_valid = size(valid_s, 1); % Number of remaining valid states

        if num_removed > 0 && ~cfg.suppress_console_output
             fprintf('  Removed %d invalid/shoot-through state(s).\n', num_removed);
        elseif num_removed == 0 && n_sw >=2 && ~cfg.suppress_console_output
             fprintf('  No invalid/shoot-through states removed (besides potentially the all-ON state if n_sw < 2).\n');
        end
        
        if num_valid == 0 && n_sw > 0
            error('AnalyzeSwitchingConverter:NoValidStatesAfterFilter', 'No valid switching states remain after filtering. Check netlist and switch definitions.');
        end

    end % end generateValidStates

    %------------------------------------------------
    function seq_out = validateStateSequence(seq_in, num_valid)
        % Validates the state sequence or sets the default.
        if isempty(seq_in)
            if num_valid == 0 % Should not happen if generateValidStates errors out
                 seq_out = [];
            else
                 seq_out = 1:num_valid; % Default sequence [1, 2, ..., num_valid]
            end
        else
            % Validate user-provided sequence
            if ~isnumeric(seq_in) || ~isvector(seq_in) || isempty(seq_in) || any(seq_in < 1) || any(seq_in > num_valid) || any(mod(seq_in,1) ~= 0)
                error('AnalyzeSwitchingConverter:InvalidStateSequence', ...
                      'config.state_sequence must be a non-empty vector of valid state indices (integers between 1 and %d).', num_valid);
            end
             seq_out = seq_in;
        end
    end % end validateStateSequence

    %% --- Phase 2 Helpers ---
    function node_map_updated = updateNodeMapForParasitics(nodeMap_in, switches_in)
        % Adds internal nodes for parasitic components (Ron/Roff || (Rs + Cs)) model.
        % Does not perform clash checks as per user confirmation.
        node_map_updated = nodeMap_in;
        map_values = values(nodeMap_in);
        if isempty(map_values)
            max_node_idx = 0;
        else
            max_node_idx = max(cell2mat(map_values)); % Find current max index
        end

        for sw_idx = 1:length(switches_in)
            sw_name = switches_in(sw_idx).name;
            % Create a unique internal node name (e.g., M1_RsCs_mid)
            internal_node_name = sprintf('%s_RsCs_mid', sw_name);

            if ~isKey(node_map_updated, internal_node_name)
                max_node_idx = max_node_idx + 1;
                node_map_updated(internal_node_name) = max_node_idx;
                fprintf('    Added parasitic node: %s (index %d)\n', internal_node_name, max_node_idx);
            else
                 % Node already exists - this case should ideally not happen with this naming
                 % but we don't error out based on user feedback.
                 warning('AnalyzeSwitchingConverter:NodeNameExists', 'Parasitic node name "%s" already exists in map. Reusing existing index.', internal_node_name);
            end
        end
    end % end updateNodeMapForParasitics

    %------------------------------------------------
    function comps_state = createStateComponentList(components_in, switches_in, state_vector, cfg, nodeMap_final_in, current_state_loop_idx)
        % Creates a temporary component list for a specific switch state.
        % Replaces 'SW' type components with 'R' (Ron/Roff).
        % If cfg.include_parasitics is true, adds RC parasitics using
        % individual values from cfg.switch_parasitics.
        % Implements Ron/Roff || (Rs + Cs) model.

        % Create a deep copy (structs are value types, cell arrays need care)
        comps_state = cell(size(components_in));
        for k_comp_copy=1:length(components_in)
            comps_state{k_comp_copy} = components_in{k_comp_copy};
        end

        switch_names = {switches_in.name};
        num_sw = length(switches_in);
        comps_to_add = {}; % Cell array for new parasitic components

        for comp_list_idx = 1:length(comps_state)
            % Use a temporary variable for the component being modified
            current_comp = comps_state{comp_list_idx};

            if isfield(current_comp,'is_switch') && current_comp.is_switch
                sw_name = current_comp.name;
                sw_idx_in_list = find(strcmp(switch_names, sw_name));

                if ~isempty(sw_idx_in_list) && sw_idx_in_list <= num_sw
                    % Get the state (ON=1, OFF=0) for this switch
                    current_switch_state = 0; % Default to OFF if no state vector (no switches case)
                    if ~isempty(state_vector)
                        current_switch_state = state_vector(sw_idx_in_list);
                    end

                    % --- Modify the original switch component to be Ron/Roff ---
                    current_comp.type = 'R';
                    current_comp.is_switch = false;
                    current_comp.is_resistor = true;
                    current_comp.is_inductor = false; current_comp.is_capacitor = false;
                    current_comp.is_voltage_source = false; current_comp.is_current_source = false;
                    current_comp.is_vcvs = false; current_comp.is_vccs = false;
                    current_comp.is_ccvs = false; current_comp.is_cccs = false;
                    current_comp.is_mutual_ind = false;

                    if current_switch_state == 1 % ON state
                        current_comp.value = cfg.Ron;
                        current_comp.name = [sw_name '_Ron']; 
                    else % OFF state
                        current_comp.value = cfg.Roff;
                        current_comp.name = [sw_name '_Roff']; 
                    end
                    current_comp.controlling_nodes = {}; % Remove control info
                    current_comp.controlling_element = '';

                    % --- Add Parasitic Components if enabled ---
                    if cfg.include_parasitics
                        % Model: Ron/Roff (already modified) in parallel with series Rs+Cs

                        n1_name = current_comp.nodes{1}; % Drain
                        n2_name = current_comp.nodes{2}; % Source
                        mid_node_name = sprintf('%s_RsCs_mid', sw_name); % Internal node

                        % Check if nodes exist in the final map (should exist if map updated correctly)
                        if ~isKey(nodeMap_final_in, n1_name) || ~isKey(nodeMap_final_in, n2_name) || ~isKey(nodeMap_final_in, mid_node_name)
                             error('AnalyzeSwitchingConverter:ParasiticNodeMissing', 'Cannot add parasitics for %s: One or more nodes (%s, %s, %s) not found in final node map.', sw_name, n1_name, n2_name, mid_node_name);
                        end

                        
                        % Retrieve switch-specific parasitic values from config struct
                        if ~isfield(cfg.switch_parasitics, sw_name)
                            % This should have been caught by validation earlier, but double-check
                            error('InternalError:ParasiticLookupFailed', 'Parasitic data for switch %s not found in config.switch_parasitics during state component list creation.', sw_name);
                        end
                        sw_parasitics = cfg.switch_parasitics.(sw_name);
                        Rs_val = sw_parasitics.Rs; % Get Rs for this specific switch
                        Cs_val = sw_parasitics.Cs; % Get Cs for this specific switch
                        

                        % 1. Rs (Series R) between Drain(N1) and MidNode
                        comp_Rs = struct();
                        comp_Rs.name = [sw_name '_Rs_parasitic']; % Unique name for parasitic R
                        comp_Rs.type = 'R';
                        comp_Rs.nodes = {n1_name, mid_node_name};
                        comp_Rs.value = Rs_val; % <<< Use switch-specific Rs
                        % Initialize all flags correctly for a resistor
                        comp_Rs.is_resistor = true; comp_Rs.is_inductor = false; comp_Rs.is_capacitor = false;
                        comp_Rs.is_voltage_source = false; comp_Rs.is_current_source = false;
                        comp_Rs.is_vcvs = false; comp_Rs.is_vccs = false; comp_Rs.is_ccvs = false;
                        comp_Rs.is_cccs = false; comp_Rs.is_mutual_ind = false; comp_Rs.is_switch = false;
                        comp_Rs.controlling_nodes = {}; comp_Rs.controlling_element = '';
                        comps_to_add{end+1} = comp_Rs; % Add to temporary list

                        % 2. Cs (Series C) between MidNode and Source(N2)
                        comp_Cs = struct();
                        comp_Cs.name = [sw_name '_Cs_parasitic']; % Unique name for parasitic C (used for ZVS target)
                        comp_Cs.type = 'C';
                        comp_Cs.nodes = {mid_node_name, n2_name};
                        comp_Cs.value = Cs_val; % <<< Use switch-specific Cs
                         % Initialize all flags correctly for a capacitor
                        comp_Cs.is_capacitor = true; comp_Cs.is_resistor = false; comp_Cs.is_inductor = false;
                        comp_Cs.is_voltage_source = false; comp_Cs.is_current_source = false;
                        comp_Cs.is_vcvs = false; comp_Cs.is_vccs = false; comp_Cs.is_ccvs = false;
                        comp_Cs.is_cccs = false; comp_Cs.is_mutual_ind = false; comp_Cs.is_switch = false;
                        comp_Cs.controlling_nodes = {}; comp_Cs.controlling_element = '';
                        comps_to_add{end+1} = comp_Cs; % Add to temporary list

                        % Original switch component (now Ron/Roff) is not removed, it stays in parallel
                    end % end if cfg.include_parasitics

                    % Update the component in the main list (it's now Ron or Roff)
                    comps_state{comp_list_idx} = current_comp;

                else
                    % This should not happen if parsing was successful
                    warning('AnalyzeSwitchingConverter:SwitchNotFound', 'Switch component %s defined in netlist not found in internal switch list.', sw_name);
                end
            end % end if is_switch
        end % end loop through original components

        % Add the new parasitic components to the end of the list for this state
        comps_state = [comps_state, comps_to_add];

        % Re-index components after additions (important for MNA)
        for k_reindex=1:length(comps_state)
            comps_state{k_reindex}.index = k_reindex;
        end

    end % end createStateComponentList


    %------------------------------------------------
    function [An, Bn, Xn, state_dim_n, mna_results] = generate_AB_for_state(components_state, nodeMap_in, coupling_map_in, cfg, u_vars_in, current_state_vector_in, current_state_idx_in)
        % Encapsulates the MNA -> A/B logic for a specific switching state.
        % Logic from netliststatespace_lcapy_Update10_4.txt.
        % Errors out if MNA matrix is degenerate.

        % --- Step 1: Identify State Variables, Inputs ---
        inductors = {}; L_values = containers.Map('KeyType', 'char', 'ValueType', 'double');
        capacitors = {}; C_values = containers.Map('KeyType', 'char', 'ValueType', 'double');
        % u_vars_in is passed as input
        num_inputs = length(u_vars_in);

        x_vars = {}; % State variables: iL, vC (Order: L first, then C)
        state_var_indices = []; % Indices in components_state

        num_L = 0; num_C = 0;

        
        % First pass: Identify L to define state variables in order
        for k_comp_L = 1:length(components_state)
            comp = components_state{k_comp_L};
            if comp.is_inductor
                num_L = num_L + 1;
                inductors{end+1} = comp.name;
                L_values(comp.name) = comp.value;
                x_vars{end+1} = ['i_' comp.name]; % Use i_ prefix for inductor current states
                state_var_indices(end+1) = k_comp_L;
            end
        end
        % Second pass: Identify C
        for k_comp_C = 1:length(components_state)
             comp = components_state{k_comp_C};
 
             if comp.is_capacitor
                num_C = num_C + 1;
                capacitors{end+1} = comp.name;
                C_values(comp.name) = comp.value;
                x_vars{end+1} = ['v_' comp.name]; % Use v_ prefix for capacitor voltage states
                state_var_indices(end+1) = k_comp_C;
            end
        end
        state_dim_n = num_L + num_C;
        Xn = x_vars; % Return the identified state variable names/order

        fprintf('      Identified %d states (%d L, %d C): %s\n', state_dim_n, num_L, num_C, strjoin(Xn, ', '));

        % Handle case with no state variables
        if state_dim_n == 0
            An = zeros(0, 0);
            Bn = zeros(0, num_inputs);
            % Create a minimal mna_results struct for derive_CD
            mna_results = struct();
            mna_results.S_states = zeros(length(keys(nodeMap_in))-1, 0); % N_mna x 0
            mna_results.S_inputs = zeros(length(keys(nodeMap_in))-1, num_inputs); % N_mna x num_inputs
            mna_results.nodeMap = nodeMap_in;
            mna_results.extra_vars_map = containers.Map(); % Empty map
            mna_results.N_mna = length(keys(nodeMap_in)) - 1; % Assume only node voltages if no V sources
            mna_results.num_nodes = mna_results.N_mna;
            return;
        end

        % Create component map for this state (mapping name to component struct)
        comp_map_state = containers.Map('KeyType', 'char', 'ValueType', 'any');
        for k_map=1:length(components_state), comp_map_state(components_state{k_map}.name) = components_state{k_map}; end

        % --- Step 2: Create Substitute Circuit ---
        sscct = create_substitute_circuit_full_k(components_state, nodeMap_in); % Nested helper

        % --- Step 3: Perform MNA on Substitute Circuit ---
        % This nested function now performs the MNA calculation
        mna_results = perform_mna_sscct_full_k(sscct, state_dim_n, num_inputs, Xn, u_vars_in, comp_map_state, cfg.degeneracy_threshold, current_state_vector_in, current_state_idx_in);

        if ~isempty(mna_results.S_states)
            fprintf('        DEBUG A/B Gen: norm(S_states) = %.3e, range = [%.3e, %.3e]\n', norm(mna_results.S_states), min(mna_results.S_states(:)), max(mna_results.S_states(:)));
        else
            fprintf('        DEBUG A/B Gen: S_states is empty.\n');
        end
        if ~isempty(mna_results.S_inputs)
             fprintf('        DEBUG A/B Gen: norm(S_inputs) = %.3e, range = [%.3e, %.3e]\n', norm(mna_results.S_inputs), min(mna_results.S_inputs(:)), max(mna_results.S_inputs(:)));
        else
             fprintf('        DEBUG A/B Gen: S_inputs is empty.\n');
        end

        % --- Step 4: Derive State Matrices A and B ---
        An = zeros(state_dim_n, state_dim_n);
        Bn = zeros(state_dim_n, num_inputs);
        processed_inductors = containers.Map('KeyType', 'char', 'ValueType', 'logical'); % Track processed coupled inductors

        % Process Inductors (including mutual coupling)
        inductor_state_indices = find(startsWith(Xn, 'i_')); % Indices of inductor current states in Xn
        for idx_in_Xn = inductor_state_indices
            state_var_name1 = Xn{idx_in_Xn};
            comp_name1 = extractAfter(state_var_name1, 'i_'); % Get original inductor name

            if isKey(processed_inductors, comp_name1) && processed_inductors(comp_name1)
                continue; % Already processed as part of a pair
            end

            comp1 = comp_map_state(comp_name1);
            L1_val = L_values(comp_name1);
            row_idx1 = idx_in_Xn; % Row in A/B matrix corresponds to index in Xn

            % Get vL1 = Vn1a - Vn1b coefficients from MNA sensitivities
            % Use the nested get_output_coeffs function
            [~, Cx_vL1, Du_vL1] = get_output_coeffs_numeric_full_k('V', comp1.name, sscct, mna_results, Xn, u_vars_in, comp_map_state);

            % Check for coupling
            coupled_inductor_name = ''; k_val = 0; M_val = 0; L2_val = 0;
            row_idx2 = []; comp2 = []; Cx_vL2 = []; Du_vL2 = []; state_var_name2 = '';
            is_coupled_pair = false;
            if isKey(coupling_map_in, comp_name1)
                mutual_info = coupling_map_in(comp_name1);
                % Determine the name of the OTHER inductor in the pair
                if strcmp(mutual_info.L1_name, comp_name1)
                    coupled_inductor_name = mutual_info.L2_name;
                else
                    coupled_inductor_name = mutual_info.L1_name;
                end
                k_val = mutual_info.k;

                % Check if the coupled inductor exists in this state's list and hasn't been processed
                if isKey(L_values, coupled_inductor_name) && (~isKey(processed_inductors, coupled_inductor_name) || ~processed_inductors(coupled_inductor_name))
                    is_coupled_pair = true;
                    comp2 = comp_map_state(coupled_inductor_name);
                    L2_val = L_values(coupled_inductor_name);
                    state_var_name2 = ['i_' coupled_inductor_name];
                    row_idx2 = find(strcmp(Xn, state_var_name2)); % Find index of coupled inductor state in Xn
                    if isempty(row_idx2)
                        error('InternalError:CoupledStateVarNotFound', 'Coupled state variable %s not found in state list Xn.', state_var_name2);
                    end
                    if L1_val <= 0 || L2_val <= 0
                        error('AnalyzeSwitchingConverter:InvalidInductanceValue', 'Self-inductance L1=%g or L2=%g is non-positive for coupled pair %s-%s.', L1_val, L2_val, comp_name1, coupled_inductor_name);
                    end
                    M_val = k_val * sqrt(L1_val * L2_val); % Calculate mutual inductance M
                    [~, Cx_vL2, Du_vL2] = get_output_coeffs_numeric_full_k('V', comp2.name, sscct, mna_results, Xn, u_vars_in, comp_map_state);
                    fprintf('      Processing coupled inductors: %s, %s (k=%.3f, M=%.3e)\n', comp_name1, coupled_inductor_name, k_val, M_val);
                end
            end

            % --- Calculate A, B entries ---
            if is_coupled_pair
                % Form L matrix [L1 M; M L2]
                L_matrix = [L1_val, M_val; M_val, L2_val];
                det_L = det(L_matrix);
                % Check singularity relative to magnitudes
                if abs(det_L) < eps(max(L1_val*L2_val, M_val^2))
                    error('AnalyzeSwitchingConverter:InvalidInductanceMatrix', 'Inductance matrix for %s-%s is singular (L1=%g, L2=%g, M=%g, det=%g). Coupling k=%.4f too close to 1?', comp_name1, coupled_inductor_name, L1_val, L2_val, M_val, det_L, k_val);
                end
                inv_L_matrix = inv(L_matrix);

                % Combine voltage coefficient vectors: [vL1_coeffs; vL2_coeffs]
                vL_coeffs_x = [Cx_vL1; Cx_vL2];
                vL_coeffs_u = [Du_vL1; Du_vL2];

                % Calculate di/dt = inv(L) * vL
                di_dt_coeffs_x = inv_L_matrix * vL_coeffs_x;
                di_dt_coeffs_u = inv_L_matrix * vL_coeffs_u;

                % Assign rows to A and B matrices
                An(row_idx1, :) = di_dt_coeffs_x(1, :);
                An(row_idx2, :) = di_dt_coeffs_x(2, :);
                if num_inputs > 0
                    Bn(row_idx1, :) = di_dt_coeffs_u(1, :);
                    Bn(row_idx2, :) = di_dt_coeffs_u(2, :);
                end

                % Mark both as processed
                processed_inductors(comp_name1) = true;
                processed_inductors(coupled_inductor_name) = true;

            else % Uncoupled Inductor
                 if L1_val <= 0, error('AnalyzeSwitchingConverter:InvalidInductanceValue', 'Invalid or non-positive inductance for %s: %g', comp_name1, L1_val); end
                % diL/dt = vL/L
                % +++ ADD THESE PRINT STATEMENTS +++
                fprintf('        DEBUG A/B Gen (L): L1=%s, L1_val=%.3e, norm(Cx_vL1)=%.3e, norm(Du_vL1)=%.3e\n', comp_name1, L1_val, norm(Cx_vL1), norm(Du_vL1));
                % +++ END ADDED PRINT STATEMENTS +++
                An(row_idx1, :) = Cx_vL1 / L1_val;

                if num_inputs > 0
                    Bn(row_idx1, :) = Du_vL1 / L1_val;
                end

                 fprintf('        DEBUG A/B Gen (L): Result norm(An(%d,:))=%.3e, norm(Bn(%d,:))=%.3e\n', row_idx1, norm(An(row_idx1, :)), row_idx1, norm(Bn(row_idx1, :)));

                processed_inductors(comp_name1) = true;
            end
        end % End inductor loop

        % Process Capacitors: dvC/dt = iC / C
        capacitor_state_indices = find(startsWith(Xn, 'v_')); % Indices of capacitor voltage states in Xn
        for idx_in_Xn = capacitor_state_indices
            state_var_name = Xn{idx_in_Xn};
            comp_name = extractAfter(state_var_name, 'v_'); % Get original capacitor name
            C_val = C_values(comp_name);
            if C_val <= 0, error('AnalyzeSwitchingConverter:InvalidCapacitanceValue', 'Invalid or non-positive capacitance for %s: %g', comp_name, C_val); end
            row_idx = idx_in_Xn; % Row in A/B matrix

            % Get iC coefficients from MNA sensitivities
            % Note: For capacitor C, we need current through the substitute V_ source in sscct
            [~, Cx_iC, Du_iC] = get_output_coeffs_numeric_full_k('I', comp_name, sscct, mna_results, Xn, u_vars_in, comp_map_state);
            fprintf('        DEBUG A/B Gen (C): C=%s, C_val=%.3e, norm(Cx_iC)=%.3e, norm(Du_iC)=%.3e\n', comp_name, C_val, norm(Cx_iC), norm(Du_iC));

            % dvC/dt = iC/C
            An(row_idx, :) = Cx_iC / C_val;
            if num_inputs > 0
                Bn(row_idx, :) = Du_iC / C_val;
            end

            fprintf('        DEBUG A/B Gen (C): Result norm(An(%d,:))=%.3e, norm(Bn(%d,:))=%.3e\n', row_idx, norm(An(row_idx, :)), row_idx, norm(Bn(row_idx, :)));


        end % End capacitor loop

      
        fprintf('      A and B matrices derived for state.\n');

    end % end generate_AB_for_state

   %------------------------------------------------
    function [Cn, Dn] = derive_CD_for_state(outputs_req, components_state, mna_results_in, state_vars_in, input_vars_in, cfg)
        % Derives C and D matrices for a given state based on requested outputs.
        % MODIFIED: Handles {'all_states'} to output all state variables.

        num_states = length(state_vars_in);
        num_inputs = length(input_vars_in);

        outputs_to_calculate = {}; % Initialize
        is_all_states_mode = false;

        % --- Check if 'all_states' mode is requested ---
        if iscell(outputs_req) && isscalar(outputs_req) && strcmp(outputs_req{1}, 'all_states')
            outputs_to_calculate = state_vars_in; % Use the state vars for this state as outputs
            is_all_states_mode = true;
            if ~isempty(outputs_to_calculate)
                fprintf('      Deriving C and D matrices for all %d state variables...\n', length(outputs_to_calculate));
            else
                 fprintf('      No state variables found for this state. Skipping C and D matrices.\n');
            end
        elseif iscell(outputs_req) && ~isempty(outputs_req)
            outputs_to_calculate = outputs_req; % Use user-defined list
             fprintf('      Deriving C and D matrices for %d specified outputs...\n', length(outputs_to_calculate));
        else
             fprintf('      No outputs requested or specified. Skipping C and D matrices.\n');
             outputs_to_calculate = {}; % Ensure it's empty
        end

        num_outputs = length(outputs_to_calculate);
        Cn = zeros(num_outputs, num_states);
        Dn = zeros(num_outputs, num_inputs);

        if num_outputs == 0
            return; % Return empty matrices if no outputs requested/needed
        end

        % Create component map for this state
        comp_map_state = containers.Map('KeyType', 'char', 'ValueType', 'any');
        for k_map_cd=1:length(components_state), comp_map_state(components_state{k_map_cd}.name) = components_state{k_map_cd}; end

        % Get substitute circuit (needed for get_output_coeffs)
        sscct = create_substitute_circuit_full_k(components_state, mna_results_in.nodeMap);

        for k_out = 1:num_outputs
           output_def = outputs_to_calculate{k_out};
           output_type = ''; target_name = '';

           % --- Parse the output definition ---
           if is_all_states_mode
               % Output definition is a state variable name like 'i_L1' or 'v_C1'
               if startsWith(output_def, 'i_')
                   output_type = 'I';
                   target_name = extractAfter(output_def, 'i_');
                   % Ensure target component exists
                   if ~isKey(comp_map_state, target_name) || ~comp_map_state(target_name).is_inductor
                        warning('AnalyzeSwitchingConverter:InvalidStateVarOutput', 'Inductor "%s" for state variable output "%s" not found. Skipping.', target_name, output_def);
                        continue;
                   end
               elseif startsWith(output_def, 'v_')
                   output_type = 'V';
                   target_name = extractAfter(output_def, 'v_');
                   % Ensure target component exists
                   if ~isKey(comp_map_state, target_name) || ~comp_map_state(target_name).is_capacitor
                        warning('AnalyzeSwitchingConverter:InvalidStateVarOutput', 'Capacitor "%s" for state variable output "%s" not found. Skipping.', target_name, output_def);
                        continue;
                   end
               else
                   warning('AnalyzeSwitchingConverter:InvalidStateVarOutput', 'Could not parse state variable "%s" as an output request. Skipping.', output_def);
                   continue;
               end
           else
               % Parse user-defined 'V(target)' or 'I(target)'
               match = regexp(output_def, '([VI])\((.+)\)', 'tokens');
               if isempty(match) || length(match{1}) ~= 2
                   warning('AnalyzeSwitchingConverter:InvalidOutputDef', 'Invalid output definition: "%s". Skipping. Use format V(nodename) or I(compname).', output_def);
                   continue;
               end
               output_type = upper(match{1}{1}); % 'V' or 'I'
               target_name = match{1}{2};       % Node name or Component name
           end

           % --- Calculate coefficients ---
           if isempty(output_type) || isempty(target_name)
                warning('AnalyzeSwitchingConverter:OutputParseError', 'Failed to parse output definition "%s". Skipping.', output_def);
                continue;
           end

           try
               % Call the adapted coefficient calculation function
               [~, Cx_k, Du_k] = get_output_coeffs_numeric_full_k(output_type, target_name, sscct, mna_results_in, state_vars_in, input_vars_in, comp_map_state); % Nested helper

               % Assign the calculated row to Cn and Dn
               if num_states > 0 && ~isempty(Cx_k); Cn(k_out, :) = Cx_k; end
               if num_inputs > 0 && ~isempty(Du_k); Dn(k_out, :) = Du_k; end

           catch ME_CD
                warning('AnalyzeSwitchingConverter:OutputError', 'Could not calculate output "%s": %s', output_def, ME_CD.message);
                % Leave the row as zeros if calculation fails
           end
        end % End loop through outputs
        if num_outputs > 0
            fprintf('      C and D matrices derived.\n');
        end

    end % end derive_CD_for_state

    %------------------------------------------------
    % --- Nested MNA and Coefficient Helpers ---
    %------------------------------------------------

    function sscct = create_substitute_circuit_full_k(components_in, nodeMap_in)
        % Creates the substitute circuit (L->I_, C->V_), copies others.
        sscct.elements = {};
        sscct.nodeMap = nodeMap_in;
        sscct.nodes = keys(nodeMap_in); % Store node names

        for k_sscct_create = 1:length(components_in)
            comp = components_in{k_sscct_create};
            new_comp = comp; % Start with original component data

            if new_comp.is_inductor
                % Replace L with Current Source I_Lname
                new_comp.type = 'I_'; % Mark as substitute current source
                new_comp.value = 0; % Value is symbolic iL(t), set to 0 for numeric MNA setup
                % Update flags
                new_comp.is_inductor = false; new_comp.is_current_source = true; new_comp.is_voltage_source = false;
                new_comp.is_resistor = false; new_comp.is_capacitor = false;
                sscct.elements{end+1} = new_comp;
            elseif new_comp.is_capacitor
                % Replace C with Voltage Source V_Cname
                new_comp.type = 'V_'; % Mark as substitute voltage source
                new_comp.value = 0; % Value is symbolic vC(t), set to 0 for numeric MNA setup
                % Update flags
                new_comp.is_capacitor = false; new_comp.is_voltage_source = true; new_comp.is_current_source = false;
                new_comp.is_resistor = false; new_comp.is_inductor = false;
                sscct.elements{end+1} = new_comp;
            elseif ~new_comp.is_mutual_ind % Copy R, V, I, E, G, H, F, SW(as R) directly
                sscct.elements{end+1} = new_comp;
            end
            % K elements are ignored in substitute circuit MNA
        end
    end % end create_substitute_circuit_full_k

    %------------------------------------------------
    function mna_results = perform_mna_sscct_full_k(sscct, num_states, num_inputs, state_vars, input_vars, comp_map_in, degeneracy_threshold, state_vector_for_debug, state_idx_for_debug)
        % Performs MNA on the substitute circuit including dependent sources.
        % Calculates sensitivities numerically: X_mna = S_states * x + S_inputs * u
        % Issues a WARNING if degeneracy detected but proceeds with calculation.

        nodeMap_local = sscct.nodeMap;
        num_nodes = length(keys(nodeMap_local)) - 1; % Exclude ground '0'

        % --- Identify MNA Unknowns ---
        % Unknowns are:
        % 1. Node voltages (excluding ground) -> indices 1 to num_nodes
        % 2. Currents through V, V_, E, H sources -> indices num_nodes+1 to N_mna
        extra_vars_map = containers.Map('KeyType', 'char', 'ValueType', 'double'); % Map: Vsrc_name -> MNA current index offset (1-based)
        num_extra_vars = 0;
        for k_mna_vars = 1:length(sscct.elements)
            comp = sscct.elements{k_mna_vars};
            % All voltage sources (V, V_, E, H) in the substitute circuit require a current unknown
            if comp.is_voltage_source
                if ~isKey(extra_vars_map, comp.name) % Avoid duplicates if source appears multiple times (shouldn't happen)
                    num_extra_vars = num_extra_vars + 1;
                    extra_vars_map(comp.name) = num_extra_vars; % Map name to index offset
                end
            end
        end
        N_mna = num_nodes + num_extra_vars; % Total number of MNA unknowns
        fprintf('        Substitute circuit MNA size: %d unknowns (%d nodes + %d V/V_/E/H currents)\n', N_mna, num_nodes, num_extra_vars);

        % --- Build MNA Matrices ---
        A_mna = zeros(N_mna, N_mna);         % MNA matrix (conductance + incidence)
        Z_states = zeros(N_mna, num_states); % RHS matrix for state variable dependencies
        Z_inputs = zeros(N_mna, num_inputs); % RHS matrix for input source dependencies

        % --- Stamp components into MNA matrices ---
        for k_stamp = 1:length(sscct.elements)
            comp = sscct.elements{k_stamp}; comp_type_mna = comp.type; comp_name = comp.name;

            % Get Node Indices (n1, n2 from MNA map)
            node_indices = zeros(1, length(comp.nodes)); valid_nodes = true;
            for n_idx = 1:length(comp.nodes)
                if isKey(nodeMap_local, comp.nodes{n_idx})
                    node_indices(n_idx) = nodeMap_local(comp.nodes{n_idx});
                else
                    warning('MNA:NodeMappingError', 'Node %s for %s not found.', comp.nodes{n_idx}, comp_name); valid_nodes = false; break;
                end
            end
            if ~valid_nodes, continue; end % Skip component if nodes not found
            n1 = node_indices(1); % MNA matrix index (1-based) or 0 for ground
            n2 = ifelse(length(node_indices)>1, node_indices(2), 0); % Handle single node comps if any

            % Get Controlling Node Indices (nc1, nc2 for E, G)
            nc1 = 0; nc2 = 0; valid_ctrl_nodes = true;
            if ~isempty(comp.controlling_nodes)
                 if isKey(nodeMap_local, comp.controlling_nodes{1}); nc1 = nodeMap_local(comp.controlling_nodes{1}); else; valid_ctrl_nodes = false; end
                 if valid_ctrl_nodes && length(comp.controlling_nodes) > 1
                     if isKey(nodeMap_local, comp.controlling_nodes{2}); nc2 = nodeMap_local(comp.controlling_nodes{2}); else; valid_ctrl_nodes = false; end
                 end
                 if ~valid_ctrl_nodes; warning('MNA:NodeMappingError', 'Ctrl node for %s not found.', comp_name); end
            end
            if ~valid_ctrl_nodes && (comp.is_vcvs || comp.is_vccs), continue; end % Skip if controlling nodes invalid

            % Get Controlling Current Index (k_ctrl for H, F)
            k_ctrl = 0; valid_ctrl_elem = true;
            if ~isempty(comp.controlling_element)
                 if isKey(extra_vars_map, comp.controlling_element) % Check if controlling element has a current variable
                     k_ctrl = num_nodes + extra_vars_map(comp.controlling_element); % Get index in MNA vector
                     % Check if controlling element is actually V, E, or H type in original map
                     if isKey(comp_map_in, comp.controlling_element)
                         ctrl_comp_orig = comp_map_in(comp.controlling_element);
                         if ~(ctrl_comp_orig.is_voltage_source || ctrl_comp_orig.is_vcvs || ctrl_comp_orig.is_ccvs)
                              warning('MNA:ControlSourceError', 'Component %s controlling %s must be V, E, or H type. Skipping component.', comp.controlling_element, comp_name);
                              valid_ctrl_elem = false;
                         end
                     else
                          warning('MNA:ControlSourceError', 'Original controlling source %s for component %s not found. Skipping component.', comp.controlling_element, comp_name);
                          valid_ctrl_elem = false;
                     end
                 else
                     warning('MNA:ControlSourceError', 'Controlling source %s for component %s not found in MNA variable map or is not a V/E/H type. Skipping component.', comp.controlling_element, comp_name);
                     valid_ctrl_elem = false;
                 end
            end
            if ~valid_ctrl_elem && (comp.is_ccvs || comp.is_cccs), continue; end % Skip if controlling element invalid

            % Get Current Variable Index k_var (for V, V_, E, H)
            k_var = 0; % Index of the MNA unknown corresponding to the current through this component
            if comp.is_voltage_source && isKey(extra_vars_map, comp_name)
                k_var = num_nodes + extra_vars_map(comp_name); % Index in MNA vector [Vn; Iextra]
            end

            % --- Apply MNA Stamps based on component type in substitute circuit ---
            input_idx = []; % Reset for each component
            switch comp_type_mna
                case 'R' % Resistor (includes Ron/Roff/Rs_parasitic)
                    if isnan(comp.value) || comp.value <= 0
                        warning('MNA:InvalidValue', 'Invalid or non-positive resistance value for %s. Skipping.', comp_name);
                        continue;
                    end
                    g = 1 / comp.value; % Conductance
                    % Add conductance to the G matrix part (upper left) of A_mna
                    if n1 > 0; A_mna(n1, n1) = A_mna(n1, n1) + g; end % Stamp G on diagonal
                    if n2 > 0; A_mna(n2, n2) = A_mna(n2, n2) + g; end % Stamp G on diagonal
                    if n1 > 0 && n2 > 0 % Stamp -G on off-diagonal
                        A_mna(n1, n2) = A_mna(n1, n2) - g;
                        A_mna(n2, n1) = A_mna(n2, n1) - g;
                    end
                case 'V' % Original independent V source
                    if k_var == 0, continue; end % Should have a current variable
                    input_idx = find(strcmp(input_vars, comp_name));
                    % Stamp B matrix part (incidence for V sources)
                    if n1 > 0; A_mna(n1, k_var) = A_mna(n1, k_var) + 1; end % KCL at node n1
                    if n2 > 0; A_mna(n2, k_var) = A_mna(n2, k_var) - 1; end % KCL at node n2
                    % Stamp C matrix part (incidence for V sources, = B')
                    if n1 > 0; A_mna(k_var, n1) = A_mna(k_var, n1) + 1; end % Voltage def row k_var
                    if n2 > 0; A_mna(k_var, n2) = A_mna(k_var, n2) - 1; end % Voltage def row k_var
                    % Stamp RHS for voltage definition: Vn1 - Vn2 = Uk (input)
                    if ~isempty(input_idx) && num_inputs > 0; Z_inputs(k_var, input_idx) = 1; end
                case 'I' % Original independent I source
                    input_idx = find(strcmp(input_vars, comp_name));
                    if ~isempty(input_idx) && num_inputs > 0
                        % Stamp RHS for KCL equations: add/subtract Ik (input)
                        if n1 > 0; Z_inputs(n1, input_idx) = Z_inputs(n1, input_idx) - 1; end % Current leaves n1
                        if n2 > 0; Z_inputs(n2, input_idx) = Z_inputs(n2, input_idx) + 1; end % Current enters n2
                    end
                case 'V_' % Substitute V source (from Capacitor C)
                     if k_var == 0, continue; end % Should have a current variable
                    state_idx = find(strcmp(state_vars, ['v_' comp_name])); % Find corresponding state variable index
                    if isempty(state_idx)
                        error('MNA:StateVarNotFound', 'State variable %s not found for substitute source %s', ['v_' comp_name], comp_name);
                    end
                    % Stamp B and C matrix parts (incidence)
                    if n1 > 0; A_mna(n1, k_var) = A_mna(n1, k_var) + 1; A_mna(k_var, n1) = A_mna(k_var, n1) + 1; end
                    if n2 > 0; A_mna(n2, k_var) = A_mna(n2, k_var) - 1; A_mna(k_var, n2) = A_mna(k_var, n2) - 1; end
                    % Stamp RHS for voltage definition: Vn1 - Vn2 = Xk (state variable)
                    if num_states > 0; Z_states(k_var, state_idx) = 1; end
                case 'I_' % Substitute I source (from Inductor L)
                    state_idx = find(strcmp(state_vars, ['i_' comp_name])); % Find corresponding state variable index
                     if isempty(state_idx)
                        error('MNA:StateVarNotFound', 'State variable %s not found for substitute source %s', ['i_' comp_name], comp_name);
                    end
                    % Stamp RHS for KCL equations: add/subtract Ik (state variable)
                    if n1 > 0 && num_states > 0; Z_states(n1, state_idx) = Z_states(n1, state_idx) - 1; end % Current leaves n1
                    if n2 > 0 && num_states > 0; Z_states(n2, state_idx) = Z_states(n2, state_idx) + 1; end % Current enters n2

                % --- Dependent Sources Stamps ---
                case 'E' % VCVS: V(n1,n2) = Gain * V(nc1,nc2)
                    if k_var == 0, continue; end % Should have a current variable
                    Gain = comp.value;
                    % Stamp B and C matrix parts (incidence)
                    if n1 > 0; A_mna(n1, k_var) = A_mna(n1, k_var) + 1; A_mna(k_var, n1) = A_mna(k_var, n1) + 1; end
                    if n2 > 0; A_mna(n2, k_var) = A_mna(n2, k_var) - 1; A_mna(k_var, n2) = A_mna(k_var, n2) - 1; end
                    % Modify voltage definition row k_var: Vn1 - Vn2 - Gain*Vnc1 + Gain*Vnc2 = 0
                    if nc1 > 0; A_mna(k_var, nc1) = A_mna(k_var, nc1) - Gain; end
                    if nc2 > 0; A_mna(k_var, nc2) = A_mna(k_var, nc2) + Gain; end
                case 'G' % VCCS: I(n1->n2) = Gm * V(nc1,nc2)
                    Gm = comp.value;
                    % Modify G matrix part (upper left) for KCL equations
                    if n1 > 0 && nc1 > 0; A_mna(n1, nc1) = A_mna(n1, nc1) + Gm; end % KCL at n1
                    if n1 > 0 && nc2 > 0; A_mna(n1, nc2) = A_mna(n1, nc2) - Gm; end % KCL at n1
                    if n2 > 0 && nc1 > 0; A_mna(n2, nc1) = A_mna(n2, nc1) - Gm; end % KCL at n2
                    if n2 > 0 && nc2 > 0; A_mna(n2, nc2) = A_mna(n2, nc2) + Gm; end % KCL at n2
                case 'H' % CCVS: V(n1,n2) = R_trans * I_ctrl
                    if k_var == 0, continue; end % Should have a current variable
                    R_trans = comp.value;
                    % Stamp B and C matrix parts (incidence)
                    if n1 > 0; A_mna(n1, k_var) = A_mna(n1, k_var) + 1; A_mna(k_var, n1) = A_mna(k_var, n1) + 1; end
                    if n2 > 0; A_mna(n2, k_var) = A_mna(n2, k_var) - 1; A_mna(k_var, n2) = A_mna(k_var, n2) - 1; end
                    % Modify voltage definition row k_var: Vn1 - Vn2 - R_trans * I_ctrl(k_ctrl) = 0
                    if k_ctrl > 0; A_mna(k_var, k_ctrl) = A_mna(k_var, k_ctrl) - R_trans; end % I_ctrl is MNA unknown k_ctrl
                case 'F' % CCCS: I(n1->n2) = Beta * I_ctrl
                    Beta = comp.value;
                    % Modify B matrix part (upper right) for KCL equations
                    if n1 > 0 && k_ctrl > 0; A_mna(n1, k_ctrl) = A_mna(n1, k_ctrl) + Beta; end % KCL at n1 depends on I_ctrl(k_ctrl)
                    if n2 > 0 && k_ctrl > 0; A_mna(n2, k_ctrl) = A_mna(n2, k_ctrl) - Beta; end % KCL at n2 depends on I_ctrl(k_ctrl)
            end % End switch comp_type_mna

        end % End loop through sscct elements

 % --- Check for degeneracy, Solve, and Debug Output ---
        local_S_states = []; % Initialize to ensure they exist for assignin and mna_results
        local_S_inputs = []; % Initialize

        if N_mna > 0 % Proceed only if the MNA system is non-empty
            mna_cond = rcond(A_mna); % Calculate reciprocal condition number
            fprintf('        MNA matrix condition number (rcond): %.3e\n', mna_cond);

            if mna_cond < degeneracy_threshold
                 warning('AnalyzeSwitchingConverter:DegenerateCircuitWarn', ...
                      'Circuit might be degenerate (ill-conditioned MNA matrix, rcond = %g < threshold %g). Results may be inaccurate. Proceeding with calculation.', mna_cond, degeneracy_threshold);
            end

            % --- Attempt to solve for sensitivities ---
            if num_states > 0
                try
                    if mna_cond < degeneracy_threshold % Use pinv if rcond is too low
                        fprintf('!!! MNA_DEBUG: Using PINV for S_states due to low rcond (%.3e) for state index %d\n', mna_cond, state_idx_for_debug);
                        local_S_states = pinv(A_mna) * Z_states;
                    else
                        local_S_states = A_mna \ Z_states;
                    end
                catch ME_solve_states
                    warning('MNA:SolveErrorStates', 'Error solving MNA system for S_states: %s. S_states will be NaN.', ME_solve_states.message);
                    local_S_states = nan(N_mna, num_states); % Fill with NaN on error
                end
            else
                local_S_states = zeros(N_mna, 0); 
            end

            if num_inputs > 0
                try
                    if mna_cond < degeneracy_threshold % Use pinv if rcond is too low
                        fprintf('!!! MNA_DEBUG: Using PINV for S_inputs due to low rcond (%.3e) for state index %d\n', mna_cond, state_idx_for_debug);
                        local_S_inputs = pinv(A_mna) * Z_inputs;
                    else
                        local_S_inputs = A_mna \ Z_inputs;
                    end
                catch ME_solve_inputs
                     warning('MNA:SolveErrorInputs', 'Error solving MNA system for S_inputs: %s. S_inputs will be NaN.', ME_solve_inputs.message);
                     local_S_inputs = nan(N_mna, num_inputs); % Fill with NaN on error
                end
            else
                local_S_inputs = zeros(N_mna, 0); 
            end
            % --- End solve attempt ---

            % --- Save MNA data to base workspace for debugging if rcond is low ---
            if mna_cond < degeneracy_threshold || mna_cond == 0 % Ensure it saves even if rcond is exactly 0
                fprintf('!!! DEBUG: Low rcond (%.3e) for state index %d (Switch Vector: [%s]). Saving MNA data to base workspace with prefix dbg_MNA_State%d_ ...\n', ...
                        mna_cond, state_idx_for_debug, num2str(state_vector_for_debug), state_idx_for_debug);
                
                prefix = sprintf('dbg_MNA_State%d_', state_idx_for_debug);
                
                assignin('base', [prefix 'A_mna'], A_mna);
                assignin('base', [prefix 'Z_states'], Z_states);
                assignin('base', [prefix 'Z_inputs'], Z_inputs);
                assignin('base', [prefix 'S_states_calculated'], local_S_states);
                assignin('base', [prefix 'S_inputs_calculated'], local_S_inputs);
                assignin('base', [prefix 'sscct_elements'], sscct.elements);
                assignin('base', [prefix 'nodeMap_MNA'], nodeMap_local);
                assignin('base', [prefix 'extra_vars_map_MNA'], extra_vars_map);
                assignin('base', [prefix 'state_vector_config'], state_vector_for_debug);
                assignin('base', [prefix 'rcond_val'], mna_cond);
                assignin('base', [prefix 'N_mna_val'], N_mna);
                assignin('base', [prefix 'num_nodes_val'], num_nodes);
                assignin('base', [prefix 'num_extra_vars_val'], num_extra_vars);
                assignin('base', [prefix 'state_vars_in_MNA_Xn'], state_vars); 
                assignin('base', [prefix 'input_vars_in_MNA_u'], input_vars); 
            end
            % --- End MNA debug data saving ---
        
        else % N_mna == 0 (empty circuit or only independent sources without nodes)
             local_S_states = zeros(0, num_states); 
             local_S_inputs = zeros(0, num_inputs);
        end

        % --- Package results ---
        mna_results.S_states = local_S_states; % Sensitivity of MNA vars to state vars
        mna_results.S_inputs = local_S_inputs; % Sensitivity of MNA vars to input vars
        mna_results.num_nodes = num_nodes; % Number of non-ground nodes
        mna_results.extra_vars_map = extra_vars_map; % Map: Vsrc_name -> MNA current index offset
        mna_results.N_mna = N_mna; % Size of MNA system
        mna_results.nodeMap = nodeMap_local; % Node map used for this MNA

    end % end perform_mna_sscct_full_k


   %------------------------------------------------
    function [Cz, Cx, Du] = get_output_coeffs_numeric_full_k(target_type, target_name, sscct_in, mna_results_in, state_vars_in, input_vars_in, comp_map_in)
        % Calculates how an output y depends on states x and inputs u using sensitivities.
        % Returns Cx and Du such that y = Cx * x + Du * u.
        % Handles outputs related to parasitic components.
        % Handles I(C) requests correctly.

        % --- Input parameters ---
        % target_type: 'V' or 'I'
        % target_name: Node name or Component name
        % sscct_in: Substitute circuit structure for the current state
        % mna_results_in: Results from perform_mna_sscct_full_k (sensitivities, maps)
        % state_vars_in: Cell array of state variable names (e.g., {'i_L1', 'v_C1'})
        % input_vars_in: Cell array of input variable names (e.g., {'V1', 'I1'})
        % comp_map_in: Map linking component names to their structs for the current state

        % --- Initialization ---
        num_states = length(state_vars_in);
        num_inputs = length(input_vars_in);
        N_mna = mna_results_in.N_mna; % Size of the MNA system (num_nodes + num_extra_vars)
        num_nodes = mna_results_in.num_nodes; % Number of non-ground nodes
        extra_vars_map_local = mna_results_in.extra_vars_map; % Map: Vsrc_name -> MNA current index offset
        nodeMap_local = mna_results_in.nodeMap; % Map: node_name -> MNA node index

        L = zeros(1, N_mna); % Selection vector L maps X_mna = [Vn; I_extra] to the desired output y
        direct_state_dependency = false; % Flag if output is directly a state variable
        direct_input_dependency = false; % Flag if output is directly an input variable
        Cx = zeros(1, num_states); % Output dependency on states (y = Cx*x + ...)
        Du = zeros(1, num_inputs); % Output dependency on inputs (y = ... + Du*u)

        % --- Determine the selection vector L based on output type and target ---
        if strcmp(target_type, 'V')
            % --- Voltage Output V(target) ---
            if isKey(nodeMap_local, target_name)
                % --- V(nodename) ---
                node_name = target_name;
                if strcmp(node_name, '0')
                    L = zeros(1, N_mna); % Voltage at ground is always 0
                else
                    node_mna_idx = nodeMap_local(node_name); % Get MNA index (1 to num_nodes)
                    if node_mna_idx > 0 && node_mna_idx <= num_nodes
                        L(node_mna_idx) = 1; % Select this node voltage from MNA solution vector
                    else
                        error('InternalError:InvalidNodeIndex', 'Invalid MNA index %d for node %s.', node_mna_idx, node_name);
                    end
                end
            elseif isKey(comp_map_in, target_name)
                % --- V(compname) ---
                comp = comp_map_in(target_name);
                if length(comp.nodes) < 2
                    error('OutputTargetError:InvalidCompNodes', 'Component %s needs at least two nodes for voltage calculation.', target_name);
                end
                n1_name = comp.nodes{1}; n2_name = comp.nodes{2};
                L_Vn1 = zeros(1, N_mna); L_Vn2 = zeros(1, N_mna);
                % Get selector for V(n1)
                if ~strcmp(n1_name, '0') && isKey(nodeMap_local, n1_name)
                    idx1 = nodeMap_local(n1_name);
                    if idx1 > 0 && idx1 <= num_nodes; L_Vn1(idx1) = 1; end
                end
                % Get selector for V(n2)
                if ~strcmp(n2_name, '0') && isKey(nodeMap_local, n2_name)
                    idx2 = nodeMap_local(n2_name);
                    if idx2 > 0 && idx2 <= num_nodes; L_Vn2(idx2) = 1; end
                end
                L = L_Vn1 - L_Vn2; % V = Vn1 - Vn2
            else
                error('OutputTargetNotFound:NodeOrCompV', 'Target "%s" for voltage output not found as node or component in the current state component list.', target_name);
            end

        elseif strcmp(target_type, 'I')
            % --- Current Output I(compname) ---
            comp_name = target_name;
            if ~isKey(comp_map_in, comp_name)
                error('OutputComponentNotFound:CompNotFoundI', 'Component "%s" requested for current output not found in the current state component list.', comp_name);
            end
            comp = comp_map_in(comp_name); % Get component struct from the current state's map

            % --- Determine current calculation method based on component type ---
            if comp.is_capacitor
                % ** CORRECTED LOGIC for I(C) **
                % Current through C is the current through the substitute V_ source in sscct.
                % The substitute source has the same name as the original capacitor.
                subst_comp_name = comp_name;
                if isKey(extra_vars_map_local, subst_comp_name)
                    k_idx_offset = extra_vars_map_local(subst_comp_name); % Get MNA current variable index offset
                    k_mna_idx = num_nodes + k_idx_offset; % Calculate full MNA index
                    if k_mna_idx > 0 && k_mna_idx <= N_mna
                        L(k_mna_idx) = 1; % Select the current i_VC from MNA solution
                    else
                        error('InternalError:InvalidCurrentIndexC', 'Invalid MNA index (%d) for current of substitute V_ source %s.', k_mna_idx, comp_name);
                    end
                else
                    % This might happen if the capacitor forms a loop with V sources or is otherwise redundant in MNA
                    warning('OutputCalculation:MissingSubstVarC', 'Could not find substitute V_ source for capacitor %s in extra vars map. Output current assumed zero.', comp_name);
                    L = zeros(1, N_mna); % Set L to zero if MNA variable not found
                end
            elseif comp.is_inductor
                % ** CORRECTED LOGIC for I(L) **
                % Current through L is a state variable (i_Lname)
                state_var_name = ['i_' comp_name];
                state_idx = find(strcmp(state_vars_in, state_var_name));
                if ~isempty(state_idx)
                    if num_states > 0; Cx(state_idx(1)) = 1; end % Output is exactly this state variable
                    direct_state_dependency = true; % Mark as direct dependency
                else
                    % This could happen if an inductor doesn't contribute to a state (e.g., series with I source)
                    warning('OutputCalculation:InductorNotState', 'Inductor %s does not appear to be a state variable. Cannot calculate current directly this way.', comp_name);
                    L = zeros(1, N_mna); % Cannot calculate current via MNA selection vector L in this case
                end
            elseif comp.is_resistor % Includes original R and Ron/Roff/Rs_parasitic
                n1 = 0; n2 = 0;
                if isKey(nodeMap_local, comp.nodes{1}); n1 = nodeMap_local(comp.nodes{1}); end
                if length(comp.nodes) > 1 && isKey(nodeMap_local, comp.nodes{2}); n2 = nodeMap_local(comp.nodes{2}); end
                if isnan(comp.value) || comp.value <= 0
                    warning('InvalidValue:ResistanceOutput', 'Invalid or non-positive resistance value for %s used in output calculation. Result will be zero.', comp.name);
                    L = zeros(1, N_mna);
                else
                    L_Vn1 = zeros(1, N_mna); if n1 > 0 && n1 <= num_nodes; L_Vn1(n1) = 1; end
                    L_Vn2 = zeros(1, N_mna); if n2 > 0 && n2 <= num_nodes; L_Vn2(n2) = 1; end
                    L = (L_Vn1 - L_Vn2) / comp.value; % I = (Vn1 - Vn2) / R
                end
            elseif comp.is_voltage_source && ~comp.is_vcvs && ~comp.is_ccvs % Original V source
                % Current is the MNA unknown associated with this voltage source
                if ~isKey(extra_vars_map_local, comp_name)
                    error('InternalError:VarMapMissingV', 'Current variable for voltage source %s not found in MNA map.', comp_name);
                end
                k_idx_offset = extra_vars_map_local(comp_name);
                k_mna_idx = num_nodes + k_idx_offset;
                if k_mna_idx > 0 && k_mna_idx <= N_mna
                    L(k_mna_idx) = 1; % Selects this current variable
                else
                    error('InternalError:InvalidCurrentIndexV', 'Invalid MNA index (%d) for current of %s.', k_mna_idx, comp_name);
                end
            elseif comp.is_current_source && ~comp.is_vccs && ~comp.is_cccs % Original I source
                % Output is the input value itself
                input_idx = find(strcmp(input_vars_in, comp_name));
                if ~isempty(input_idx)
                    if num_inputs > 0; Du(input_idx(1)) = 1; end % Output is exactly this input
                    direct_input_dependency = true; % Mark as direct dependency
                else
                    % This case should not happen if I sources are correctly added to u_vars
                    error('InternalError:InputVarNotFound', 'Input variable %s not found for I source.', comp_name);
                end
            elseif comp.is_vcvs % E source (VCVS)
                 % Current is the MNA unknown associated with this source
                 if ~isKey(extra_vars_map_local, comp_name)
                    error('InternalError:VarMapMissingE', 'Current variable for VCVS %s not found in MNA map.', comp_name);
                 end
                 k_idx_offset = extra_vars_map_local(comp_name);
                 k_mna_idx = num_nodes + k_idx_offset;
                 if k_mna_idx > 0 && k_mna_idx <= N_mna
                    L(k_mna_idx) = 1; % Select the current i_E
                 else
                    error('InternalError:InvalidCurrentIndexE', 'Invalid MNA index (%d) for current of %s.', k_mna_idx, comp_name);
                 end
            elseif comp.is_vccs % G source (VCCS)
                 % Current is calculated from controlling voltage: I = Gm * (Vnc1 - Vnc2)
                 Gm = comp.value;
                 nc1_name = comp.controlling_nodes{1}; nc1_idx = 0;
                 nc2_name = comp.controlling_nodes{2}; nc2_idx = 0;
                 if isKey(nodeMap_local, nc1_name); nc1_idx = nodeMap_local(nc1_name); end
                 if isKey(nodeMap_local, nc2_name); nc2_idx = nodeMap_local(nc2_name); end
                 L_Vnc1 = zeros(1, N_mna); if nc1_idx > 0 && nc1_idx <= num_nodes; L_Vnc1(nc1_idx) = 1; end
                 L_Vnc2 = zeros(1, N_mna); if nc2_idx > 0 && nc2_idx <= num_nodes; L_Vnc2(nc2_idx) = 1; end
                 L = (L_Vnc1 - L_Vnc2) * Gm;
            elseif comp.is_ccvs % H source (CCVS)
                 % Current is the MNA unknown associated with this source
                 if ~isKey(extra_vars_map_local, comp_name)
                    error('InternalError:VarMapMissingH', 'Current variable for CCVS %s not found in MNA map.', comp_name);
                 end
                 k_idx_offset = extra_vars_map_local(comp_name);
                 k_mna_idx = num_nodes + k_idx_offset;
                 if k_mna_idx > 0 && k_mna_idx <= N_mna
                    L(k_mna_idx) = 1; % Select the current i_H
                 else
                    error('InternalError:InvalidCurrentIndexH', 'Invalid MNA index (%d) for current of %s.', k_mna_idx, comp_name);
                 end
            elseif comp.is_cccs % F source (CCCS)
                 % Current is calculated from controlling current: I = Beta * I_ctrl
                 Beta = comp.value;
                 ctrl_elem_name = comp.controlling_element; % Name of the V/E/H source providing I_ctrl
                 if ~isKey(extra_vars_map_local, ctrl_elem_name)
                      error('InternalError:VarMapMissingFctrl', 'Controlling current variable for %s (needed by %s) not found in MNA map.', ctrl_elem_name, comp_name);
                 end
                 k_ctrl_offset = extra_vars_map_local(ctrl_elem_name);
                 k_ctrl_mna_idx = num_nodes + k_ctrl_offset; % Index of I_ctrl in MNA vector
                 if k_ctrl_mna_idx > 0 && k_ctrl_mna_idx <= N_mna
                     L(k_ctrl_mna_idx) = Beta; % Select I_ctrl and multiply by Beta
                 else
                     error('InternalError:InvalidCurrentIndexFctrl', 'Invalid MNA index (%d) for controlling current of %s.', k_ctrl_mna_idx, ctrl_elem_name);
                 end
            else
                % This case should not be reached if all component types are handled
                error('OutputTypeNotSupported:CurrentCalc', 'Cannot calculate current for component type of %s', comp_name);
            end
        else
            error('OutputTypeNotSupported:UnknownType', 'Unknown output type requested: %s', target_type);
        end

        % --- Calculate Cx and Du using sensitivities if not direct dependency ---
        if direct_state_dependency || direct_input_dependency
            Cz = zeros(1, N_mna); % Not needed but return for consistency if ever used
            return; % Cx or Du already set
        end

        % Otherwise, relate the MNA selection vector L to states (x) and inputs (u)
        % y = L * X_mna = L * (S_states * x + S_inputs * u)
        % Cx = L * S_states
        % Du = L * S_inputs
        if N_mna > 0 % Proceed only if MNA system exists
            S_states = mna_results_in.S_states;
            S_inputs = mna_results_in.S_inputs;

            % Calculate Cx
            if ~isempty(S_states) && size(L, 2) == size(S_states, 1)
                if num_states > 0; Cx = L * S_states; end
            elseif num_states == 0 % Handle case with no states
                Cx = zeros(1, 0);
            else % Dimension mismatch
                 if ~isempty(S_states) % Avoid error if S_states is empty due to N_mna=0
                    error('InternalError:DimensionMismatchCx', 'Dimension mismatch for Cx calculation: size(L, 2)=%d, size(S_states, 1)=%d', size(L,2), size(S_states,1));
                 end
            end

            % Calculate Du
            if ~isempty(S_inputs) && size(L, 2) == size(S_inputs, 1)
                if num_inputs > 0; Du = L * S_inputs; end
            elseif num_inputs == 0 % Handle case with no inputs
                Du = zeros(1, 0);
            else % Dimension mismatch
                 if ~isempty(S_inputs) % Avoid error if S_inputs is empty due to N_mna=0
                    error('InternalError:DimensionMismatchDu', 'Dimension mismatch for Du calculation: size(L, 2)=%d, size(S_inputs, 1)=%d', size(L,2), size(S_inputs,1));
                 end
            end
        else % N_mna == 0, circuit was empty or only independent sources?
            Cx = zeros(1, num_states); % Output must be zero if no circuit elements
            Du = zeros(1, num_inputs);
        end

        Cz = L; % Return the MNA selection vector itself (mostly for debugging)

    end % end get_output_coeffs_numeric_full_k

    %------------------------------------------------
        function value = parse_value_simple(valueStr)
        % Parses numeric values with optional standard SI suffixes, ignoring trailing units.
        % Handles integers, decimals, scientific notation (e.g., 1.2e-3).
        % Case-insensitive for suffixes and units.

        valueStr = strtrim(valueStr); % Remove leading/trailing whitespace
        value = NaN; % Default to NaN

        % Regex to capture:
        % 1: Optional sign (+ or -)
        % 2: Numeric part (digits, optional decimal, more digits)
        % 3: Optional exponent part (e.g., e-3, E+6)
        % 4: Optional SI suffix (f, p, n, u, m, k, meg, g, t) - case insensitive
        % 5: Optional trailing non-numeric characters (units like V, Hz, Ohm etc.) - these are ignored
        pattern = '^([+\-]?)(\d+\.?\d*|\.\d+)([eE][+\-]?\d+)?([fpnumkKMGT]|meg)?([a-zA-Z]*)?$';
        tokens = regexp(valueStr, pattern, 'tokens', 'ignorecase');

        if ~isempty(tokens)
            % Combine sign, number, exponent parts
            numPartStr = [tokens{1}{1}, tokens{1}{2}, tokens{1}{3}];
            siSuffix = tokens{1}{4};
            % trailingUnits = tokens{1}{5}; % Captured but ignored

            numVal = str2double(numPartStr);

            if ~isnan(numVal) && isreal(numVal)
                multiplier = 1.0;
                if ~isempty(siSuffix)
                    % Determine multiplier based on SI suffix
                    switch lower(siSuffix)
                        case 'f', multiplier = 1e-15;
                        case 'p', multiplier = 1e-12;
                        case 'n', multiplier = 1e-9;
                        case 'u', multiplier = 1e-6;
                        case 'm', multiplier = 1e-3;
                        case 'k', multiplier = 1e3;
                        case 'meg', multiplier = 1e6; % Handle 'meg' specifically
                        case 'g', multiplier = 1e9;
                        case 't', multiplier = 1e12;
                         % Note: 'M' is handled by the 'm' case in lower() if not 'meg'
                         % Check for standalone M if meg wasn't matched
                        % case 'M' % Check explicitly if not captured by 'meg' or 'm'
                        %      if ~strcmpi(siSuffix,'meg') && ~strcmpi(siSuffix,'m')
                        %          multiplier = 1e6;
                        %      end
                    end
                end
                value = numVal * multiplier; % Calculate final value
            end
        end

        % Fallback: If regex fails (e.g., plain number), try simple str2double
        if isnan(value)
            num_direct = str2double(valueStr);
             if ~isnan(num_direct) && isreal(num_direct)
                value = num_direct;
             end
        end

         % Optional: Add a warning if parsing still failed
         if isnan(value)
             % fprintf('DEBUG: parse_value_simple_robust failed for: "%s"\n', valueStr);
             warning('AnalyzeSwitchingConverter:ValueParseFail', 'Could not parse value: "%s"', valueStr);
         end
    end % end parse_value_simple_robust



    %% --- Phase 3 Helpers ---
    function topoInfo = createTopologyInfo(valid_s, sws_parsed, state_vars_final, prms)
        % createTopologyInfo: Creates topology information.
        % MODIFIED:
        % - Expects user-defined complementary switch pairs via prms.user_complementary_switch_pairs.
        % - Populates 'complementary_switch_map': Map from switch index to its complementary switch index.
        % - Removes automatic leg detection.
        % - Retains 'zvs_target_map_switch_to_cap_idx'.
        % - Retains identification of 'deadtime_states' (global all-OFF).
        % - Retains 'state_sequence' and 'deadtime_end_seq_indices' based on prms.state_sequence.
        % - Retains the old 'zvs_target_map'.

        if ~prms.suppress_console_output
            fprintf('  Creating topology information...\n');
        end
        topoInfo = struct();
        topoInfo.num_switches = prms.num_switches; 
        topoInfo.num_states = prms.num_valid_states; 
        topoInfo.state_dim = prms.state_dim;
        
        topoInfo.state_sequence = prms.state_sequence; 
        if ~prms.suppress_console_output && ~isempty(topoInfo.state_sequence)
            fprintf('    User-defined state_sequence being used for legacy ZVS mapping: %s\n', mat2str(topoInfo.state_sequence));
        elseif isempty(topoInfo.state_sequence) && ~prms.suppress_console_output
             fprintf('    User-defined state_sequence (prms.state_sequence) is empty.\n');
        end

        % --- 1. Find all-OFF state indices (used as global deadtime states) ---
        all_off_indices = [];
        if ~isempty(valid_s) 
            all_off_indices = find(all(valid_s == 0, 2)); 
        end
        topoInfo.deadtime_states = all_off_indices(:)'; 
        if ~prms.suppress_console_output
            if isempty(topoInfo.deadtime_states)
                fprintf('    No all-OFF deadtime states found in valid_states matrix.\n');
            else
                fprintf('    Identified all-OFF deadtime state index(es): %s\n', mat2str(topoInfo.deadtime_states));
            end
        end

        % --- 2. Create zvs_target_map_switch_to_cap_idx ---
        zvs_target_map_sw_to_cap = containers.Map('KeyType', 'double', 'ValueType', 'double');
        if prms.include_parasitics && prms.num_switches > 0 && ~isempty(state_vars_final)
            if ~isfield(prms, 'switches_ordered_list') || length(prms.switches_ordered_list) ~= prms.num_switches
                error('CreateTopologyInfo:SwitchesOrderedListMismatch', ...
                      'prms.switches_ordered_list is missing or its length does not match prms.num_switches.');
            end
            for i_sw = 1:prms.num_switches
                switch_name = prms.switches_ordered_list(i_sw).name;
                target_cap_name = [switch_name '_Cs_parasitic'];
                target_state_var_name = ['v_' target_cap_name];
                
                found_cap_var_idx = find(strcmp(state_vars_final, target_state_var_name));
                
                if ~isempty(found_cap_var_idx)
                    zvs_target_map_sw_to_cap(i_sw) = found_cap_var_idx(1);
                    if ~prms.suppress_console_output
                        fprintf('    Mapping switch %s (index %d) to ZVS capacitor state var: %s (index %d)\n', ...
                                switch_name, i_sw, target_state_var_name, found_cap_var_idx(1));
                    end
                else
                    zvs_target_map_sw_to_cap(i_sw) = NaN; 
                    if ~prms.suppress_console_output
                         fprintf('    No ZVS capacitor state var found for switch %s (index %d) (expected: %s).\n', ...
                                switch_name, i_sw, target_state_var_name);
                    end
                end
            end
        else
            if ~prms.suppress_console_output
                fprintf('    Parasitics not included or no switches/state_vars; zvs_target_map_switch_to_cap_idx will be empty or all NaNs.\n');
            end
            for i_sw_nan = 1:prms.num_switches
                 zvs_target_map_sw_to_cap(i_sw_nan) = NaN;
            end
        end
        topoInfo.zvs_target_map_switch_to_cap_idx = zvs_target_map_sw_to_cap;

        % --- 3. Identify Complementary Switch Pairs (Legs) from User Input ---
        complementary_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
        num_sw_local = prms.num_switches;
        switches_info_local = prms.switches_ordered_list;

        % Initialize map with NaNs
        for i_sw_init = 1:num_sw_local
            complementary_map(i_sw_init) = NaN;
        end

        if isfield(prms, 'user_complementary_switch_pairs') && ~isempty(prms.user_complementary_switch_pairs)
            if ~iscell(prms.user_complementary_switch_pairs)
                error('CreateTopologyInfo:InvalidUserPairsType', 'prms.user_complementary_switch_pairs must be a cell array of pairs.');
            end
            if ~prms.suppress_console_output
                fprintf('    Processing user-defined complementary switch pairs...\n');
            end

            % Create a quick lookup map from switch name to index in switches_ordered_list
            switch_name_to_idx_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
            for k_sw_map = 1:num_sw_local
                switch_name_to_idx_map(switches_info_local(k_sw_map).name) = k_sw_map;
            end

            for k_pair = 1:length(prms.user_complementary_switch_pairs)
                current_pair_cell = prms.user_complementary_switch_pairs{k_pair};
                if ~iscell(current_pair_cell) || length(current_pair_cell) ~= 2
                    error('CreateTopologyInfo:InvalidPairFormat', 'Pair %d in user_complementary_switch_pairs is not a cell array of two switch names.', k_pair);
                end
                sw_name1 = current_pair_cell{1};
                sw_name2 = current_pair_cell{2};

                if ~ischar(sw_name1) || ~ischar(sw_name2)
                    error('CreateTopologyInfo:InvalidSwitchNameInPair', 'Switch names in pair %d must be character strings.', k_pair);
                end
                if strcmp(sw_name1, sw_name2)
                    error('CreateTopologyInfo:SwitchPairedWithItself', 'Switch %s in pair %d cannot be paired with itself.', sw_name1, k_pair);
                end

                if ~isKey(switch_name_to_idx_map, sw_name1)
                    error('CreateTopologyInfo:SwitchName1NotFound', 'Switch name "%s" from user-defined pair %d not found in the netlist switches.', sw_name1, k_pair);
                end
                if ~isKey(switch_name_to_idx_map, sw_name2)
                    error('CreateTopologyInfo:SwitchName2NotFound', 'Switch name "%s" from user-defined pair %d not found in the netlist switches.', sw_name2, k_pair);
                end

                idx1 = switch_name_to_idx_map(sw_name1);
                idx2 = switch_name_to_idx_map(sw_name2);

                % Check if already mapped (consistency check)
                if (isKey(complementary_map, idx1) && ~isnan(complementary_map(idx1)) && complementary_map(idx1) ~= idx2) || ...
                   (isKey(complementary_map, idx2) && ~isnan(complementary_map(idx2)) && complementary_map(idx2) ~= idx1)
                    error('CreateTopologyInfo:InconsistentPairing', 'Switch %s (idx %d) or %s (idx %d) from pair %d is already part of a different complementary pair.', sw_name1, idx1, sw_name2, idx2, k_pair);
                end
                
                complementary_map(idx1) = idx2;
                complementary_map(idx2) = idx1;
                if ~prms.suppress_console_output
                    fprintf('      User-defined leg: Switch %s (idx %d) and Switch %s (idx %d).\n', sw_name1, idx1, sw_name2, idx2);
                end
            end
        else
            if prms.run_optimization && num_sw_local > 0 && ~prms.suppress_console_output
                warning('CreateTopologyInfo:NoUserPairsForZVS', 'ZVS optimization is ON, but prms.user_complementary_switch_pairs was not provided. Leg-specific deadtime logic for ZVS might not function correctly without this information.');
            elseif ~prms.suppress_console_output
                 fprintf('    User-defined complementary switch pairs not provided. complementary_switch_map will be NaNs.\n');
            end
        end
        topoInfo.complementary_switch_map = complementary_map;


        % --- 4. Identify deadtime occurrences in the USER-DEFINED state_sequence (prms.state_sequence) ---
        deadtime_end_seq_indices_out = []; 
        num_deadtimes_in_user_sequence = 0; 

        if ~isempty(topoInfo.deadtime_states) && ~isempty(prms.state_sequence)
             num_seq_steps_local = length(prms.state_sequence); 
             for k_seq_local = 1:num_seq_steps_local 
                 current_state_idx_local = prms.state_sequence(k_seq_local);
                 if ismember(current_state_idx_local, topoInfo.deadtime_states)
                     num_deadtimes_in_user_sequence = num_deadtimes_in_user_sequence + 1;
                     deadtime_end_seq_indices_out(end+1) = k_seq_local; 
                 end
             end
        end
        topoInfo.deadtime_end_seq_indices = deadtime_end_seq_indices_out; 
        if ~prms.suppress_console_output
            if num_deadtimes_in_user_sequence > 0
                fprintf('    Found %d deadtime occurrences in the user-defined state_sequence.\n', num_deadtimes_in_user_sequence);
            else
                fprintf('    Found 0 deadtime occurrences in the user-defined state_sequence.\n');
            end
        end

        % --- 5. Create OLD zvs_target_map (occurrence in user_state_sequence -> cap_state_idx) ---
        old_zvs_target_map = containers.Map('KeyType', 'double', 'ValueType', 'double'); 
        if prms.run_optimization && num_deadtimes_in_user_sequence > 0 && ~isempty(prms.state_sequence) && num_sw_local > 0
             switches_to_target_by_name = prms.zvs_switches_to_target; 
             if ~prms.suppress_console_output
                fprintf('    Analyzing user-defined state_sequence for OLD ZVS target mapping (Targeting switches by name: %s)...\n', strjoin(switches_to_target_by_name, ', '));
             end
             
             num_seq_steps_zvs_old = length(prms.state_sequence);
             for i_dt_occurrence = 1:num_deadtimes_in_user_sequence
                 k_seq_zvs_old = topoInfo.deadtime_end_seq_indices(i_dt_occurrence); 
                 current_state_idx_zvs_old = prms.state_sequence(k_seq_zvs_old);
                 next_seq_step_idx_old = mod(k_seq_zvs_old, num_seq_steps_zvs_old) + 1;
                 next_state_idx_old = prms.state_sequence(next_seq_step_idx_old);

                 if next_state_idx_old < 1 || next_state_idx_old > size(valid_s, 1) || ...
                    current_state_idx_zvs_old < 1 || current_state_idx_zvs_old > size(valid_s, 1)
                     warning('CreateTopologyInfo:InvalidStateIdxZVS_OldMap', 'Invalid state index during OLD ZVS mapping for DT occurrence %d.', i_dt_occurrence);
                     old_zvs_target_map(i_dt_occurrence) = NaN;
                     continue;
                 end
                 
                 current_state_vector_old = valid_s(current_state_idx_zvs_old, :);
                 next_state_vector_old = valid_s(next_state_idx_old, :);
                 turning_on_switch_indices_in_valid_s = find(current_state_vector_old == 0 & next_state_vector_old == 1);
                 
                 target_cap_state_idx_for_occurrence = NaN;
                 if ~isempty(turning_on_switch_indices_in_valid_s)
                     for k_on_sw_col = 1:length(turning_on_switch_indices_in_valid_s)
                         actual_col_idx_in_valid_s = turning_on_switch_indices_in_valid_s(k_on_sw_col);
                         if actual_col_idx_in_valid_s <= length(prms.switches_ordered_list)
                             turning_on_switch_name = prms.switches_ordered_list(actual_col_idx_in_valid_s).name;
                             if ismember(turning_on_switch_name, switches_to_target_by_name)
                                 if isKey(zvs_target_map_sw_to_cap, actual_col_idx_in_valid_s)
                                     cap_var_idx_of_turning_on_sw = zvs_target_map_sw_to_cap(actual_col_idx_in_valid_s);
                                     if ~isnan(cap_var_idx_of_turning_on_sw)
                                         target_cap_state_idx_for_occurrence = cap_var_idx_of_turning_on_sw;
                                         if ~prms.suppress_console_output
                                             fprintf('        OLD MAP: DT Occurrence %d (user_seq pos %d) -> Switch %s (index %d) turning ON -> Cap Var Idx %d.\n', ...
                                                     i_dt_occurrence, k_seq_zvs_old, turning_on_switch_name, actual_col_idx_in_valid_s, target_cap_state_idx_for_occurrence);
                                         end
                                         break; 
                                     end
                                 end
                             end
                         end
                     end
                 end
                 old_zvs_target_map(i_dt_occurrence) = target_cap_state_idx_for_occurrence;
             end
        end
        topoInfo.zvs_target_map = old_zvs_target_map; 

        if ~prms.suppress_console_output
            fprintf('  Topology information created.\n');
        end
    end % end createTopologyInfo

    %------------------------------------------------
    function params_out = initializeDeadtimes(params_in)
        % initializeDeadtimes: Sets up the .dt field in params_out with per-switch deadtimes.
        % - Expects params_in.initial_deadtimes to be a vector of per-switch deadtimes.
        % - Validates its length against params_in.num_switches.
        % - Directly assigns params_in.initial_deadtimes to params_out.dt.
        % - Old logic for sequence-based deadtime interval counting is removed.
        % - Does not disable optimization here.

        params_out = params_in; % Start with input parameters
        
        if ~params_out.suppress_console_output
            fprintf('    Initializing .dt with per-switch deadtimes from .initial_deadtimes field.\n');
        end

        num_sw = params_in.num_switches;
        initial_dt_per_switch = params_in.initial_deadtimes; % This is config.initial_deadtimes from runner

        if num_sw == 0
            if ~isempty(initial_dt_per_switch)
                warning('InitializeDeadtimes:DeadtimesWithNoSwitches', ...
                        'Initial deadtimes provided, but no switches found. Setting params_out.dt to empty.');
                params_out.dt = [];
            else
                params_out.dt = []; % Correctly empty if no switches and no initial deadtimes
                if ~params_out.suppress_console_output
                    fprintf('    No switches found. params_out.dt set to empty.\n');
                end
            end
        else % num_sw > 0
            if length(initial_dt_per_switch) == num_sw
                params_out.dt = initial_dt_per_switch(:)'; % Ensure row vector
                if ~params_out.suppress_console_output
                    fprintf('    Using initial per-switch deadtimes for %d switches: [%s] ns\n', ...
                            num_sw, sprintf('%.3f ', params_out.dt*1e9));
                end
            else
                error('InitializeDeadtimes:PerSwitchDtLengthMismatch', ...
                      'Length of initial_deadtimes (%d) must match number of switches (%d).', ...
                      length(initial_dt_per_switch), num_sw);
            end
        end

        % run_optimization flag is NOT changed here. findOptimalDeadtimes will decide if it can run.
        
        if ~params_out.suppress_console_output
            if params_out.run_optimization
                fprintf('    Optimization is ON. The .dt values will be used as initial guesses by findOptimalDeadtimes.\n');
            else
                fprintf('    Optimization is OFF. The .dt values will be used as fixed applied deadtimes.\n');
            end
            fprintf('    params_out.dt initialized to: [%s]\n', mat2str(params_out.dt));
        end

    end % end initializeDeadtimes


    %% --- Phase 4 Helpers (Optimization) ---

  function state_vectors_per_step = simulateSequenceBehavior(prms, x0, As, Bs, intervals)
        % Simulates one switching cycle step-by-step.
        % Uses durations from 'intervals' for non-deadtime states
        % and specific values from 'prms.dt' for deadtime states.
        % Returns state vectors at the end of EACH STEP in the sequence.

        state_sequence = prms.state_sequence;
        num_seq_steps = length(state_sequence);
        state_dim_local = prms.state_dim;
        u_vec = constructInputVector(prms); % Get input vector
        deadtime_states_indices = prms.topology_info.deadtime_states;
        dt_vector = prms.dt; % The vector of current deadtime values
        dt_vector_idx = 1; % Index to track which deadtime value to use

        state_vectors_per_step = cell(1, num_seq_steps);
        current_state_vec = x0; % Start with initial condition

        for k_seq = 1:num_seq_steps
            state_idx = state_sequence(k_seq); % State index for this step

            % Check if state index and models are valid
            if state_idx < 1 || state_idx > prms.num_valid_states || ...
               state_idx > length(As) || isempty(As{state_idx}) || ...
               state_idx > length(Bs) || isempty(Bs{state_idx})
                 error('SimSeqBehavior:InvalidStateIdx', 'Invalid state index %d or missing model at sequence step %d.', state_idx, k_seq);
            end

            A = As{state_idx};
            B = Bs{state_idx};

            % --- Determine correct duration for this step ---
            dt_interval = 0; % Default duration
            is_deadtime_step = ismember(state_idx, deadtime_states_indices);

            if is_deadtime_step
                % Use the specific deadtime value from dt_vector
                if dt_vector_idx <= length(dt_vector)
                    dt_interval = dt_vector(dt_vector_idx);
                    dt_vector_idx = dt_vector_idx + 1; % Increment for next deadtime
                else
                    warning('SimSeqBehavior:DtVectorMismatch', 'More deadtime states in sequence than values in prms.dt at step %d.', k_seq);
                    dt_interval = 0; % Use zero if dt vector runs out
                end
            else
                % Use the pre-calculated base interval for non-deadtime states
                 if state_idx <= length(intervals)
                    dt_interval = intervals(state_idx);
                 else
                      error('SimSeqBehavior:InvalidIntervalIdx', 'Invalid state index %d for intervals array at sequence step %d.', state_idx, k_seq);
                 end
            end
            % --- End Duration Determination ---

            % Propagate state only if duration is significant
            if dt_interval > 1e-20 % Use a small threshold
                lambda_k = calculateLambda(A, dt_interval, prms); % Use existing helper
                Phi_k = eye(state_dim_local) + A * lambda_k;
                if size(B, 2) == prms.num_inputs
                    Gamma_k = lambda_k * B * u_vec;
                else
                     error('SimSeqBehavior:DimensionMismatchB', 'Dimension mismatch B matrix state %d.', state_idx);
                end
                current_state_vec = Phi_k * current_state_vec + Gamma_k; % Update state
            % else: If interval is zero, state doesn't change, current_state_vec remains the same
            end

            % Store state vector at the end of this sequence step
            state_vectors_per_step{k_seq} = current_state_vec;
        end
    end % end simulateSequenceBehavior

      %------------------------------------------------------------------
function [capacitor_voltages_per_switch, voltage_rates_per_switch] = simulateZVSConditions(prms_zvs, x0_iter, As, Bs)
        % simulateZVSConditions: Simulates one switching cycle based on dynamically
        % generated sequence and determines target capacitor voltages and rates
        % for ZVS optimization, on a per-switch basis.
        % Uses complementary_switch_map to identify leg-specific deadtimes.
        % Falls back to global deadtime_states if no complement is defined for a switch.
        % Penalty voltage in error condition is now based on primary input voltage.

        if ~prms_zvs.suppress_console_output
            fprintf('    simulateZVSConditions: Evaluating ZVS conditions for current deadtimes: [%s] ns\n', sprintf('%.3f ', prms_zvs.deadtimes_applied*1e9));
        end

        % Initialize output arrays 
        capacitor_voltages_per_switch = NaN(1, prms_zvs.num_switches);
        voltage_rates_per_switch = NaN(1, prms_zvs.num_switches);

        if prms_zvs.state_dim == 0
            if ~prms_zvs.suppress_console_output; fprintf('    simulateZVSConditions: No state variables. Skipping ZVS condition simulation.\n'); end
            return; 
        end
        
        if ~any(prms_zvs.is_switch_targeted_for_zvs)
             if ~prms_zvs.suppress_console_output; fprintf('    simulateZVSConditions: No switches are targeted for ZVS in this configuration.\n'); end
            return; 
        end

        % --- Step 1: Generate the dynamic sequence ---
        prms_for_seq_gen_zvs = prms_zvs;
        if ~isfield(prms_zvs, 'deadtimes_applied') && isfield(prms_zvs, 'dt') 
            prms_for_seq_gen_zvs.deadtimes_applied = prms_zvs.dt;
        elseif ~isfield(prms_for_seq_gen_zvs, 'deadtimes_applied')
             error('simulateZVSConditions:MissingDeadtimesApplied', 'prms_zvs.deadtimes_applied is missing for sequence generation.');
        end

        [actual_state_sequence_zvs, actual_interval_durations_zvs] = generateSequenceFromTimings(prms_for_seq_gen_zvs);

        if isempty(actual_state_sequence_zvs)
            if ~prms_zvs.suppress_console_output
                warning('simulateZVSConditions:EmptyDynamicSequence', 'generateSequenceFromTimings returned an empty sequence. Assigning penalties to targeted switches.');
            end
            for j_sw_empty_seq = 1:prms_zvs.num_switches
                if prms_zvs.is_switch_targeted_for_zvs(j_sw_empty_seq) && ~isnan(prms_zvs.zvs_cap_state_var_indices(j_sw_empty_seq))
                    penalty_V = 100.0; % Default large voltage
                    if ~isempty(prms_zvs.u_vars) && ~isempty(prms_zvs.input_values) && isKey(prms_zvs.input_values, prms_zvs.u_vars{1})
                        val = prms_zvs.input_values(prms_zvs.u_vars{1});
                        if isscalar(val) && isnumeric(val) && abs(val) > 1e-9 
                            penalty_V = abs(val);
                        end
                    end
                    capacitor_voltages_per_switch(j_sw_empty_seq) = penalty_V;
                    voltage_rates_per_switch(j_sw_empty_seq) = penalty_V / (0.05 * prms_zvs.T_period); 
                     if ~prms_zvs.suppress_console_output
                        fprintf('    simulateZVSConditions: PENALTY (Empty Sequence) for Switch %d (%s): Vc=%.4e, dVcdt=%.4e\n', ...
                            j_sw_empty_seq, prms_zvs.switches_ordered_list(j_sw_empty_seq).name, ...
                            capacitor_voltages_per_switch(j_sw_empty_seq), voltage_rates_per_switch(j_sw_empty_seq));
                    end
                end
            end
            return; 
        end
        num_dynamic_steps = length(actual_state_sequence_zvs);

        % --- Step 2: Simulate one full cycle of this dynamic sequence ---
        state_vectors_at_interval_ends = simulate_dynamic_sequence_step_ends(prms_zvs, x0_iter, As, Bs, actual_state_sequence_zvs, actual_interval_durations_zvs);
        
        u_vec_zvs = constructInputVector(prms_zvs); 
        complementary_map_local = prms_zvs.topology_info.complementary_switch_map; % Get the map

        % --- Step 3: Evaluate ZVS conditions for each targeted switch ---
        for j_sw = 1:prms_zvs.num_switches
            if prms_zvs.is_switch_targeted_for_zvs(j_sw)
                target_cap_var_idx = prms_zvs.zvs_cap_state_var_indices(j_sw);
                
                if isnan(target_cap_var_idx)
                    if ~prms_zvs.suppress_console_output
                        fprintf('    simulateZVSConditions: ZVS Cap Idx is NaN for targeted Switch %d (%s). Assigning PENALTY.\n', j_sw, prms_zvs.switches_ordered_list(j_sw).name);
                    end
                    penalty_V_nan_idx = 100.0; 
                     if ~isempty(prms_zvs.u_vars) && ~isempty(prms_zvs.input_values) && isKey(prms_zvs.input_values, prms_zvs.u_vars{1})
                        val_nan_idx = prms_zvs.input_values(prms_zvs.u_vars{1});
                        if isscalar(val_nan_idx) && isnumeric(val_nan_idx) && abs(val_nan_idx) > 1e-9
                            penalty_V_nan_idx = abs(val_nan_idx);
                        end
                    end
                    capacitor_voltages_per_switch(j_sw) = penalty_V_nan_idx; 
                    voltage_rates_per_switch(j_sw) = penalty_V_nan_idx / (0.05 * prms_zvs.T_period); 
                    continue; 
                end

                found_zvs_interval_for_switch = false;
                for k_interval = 1:num_dynamic_steps 
                    current_interval_state_idx = actual_state_sequence_zvs(k_interval);
                    current_interval_switch_config = prms_zvs.valid_states(current_interval_state_idx, :);

                    if (k_interval < num_dynamic_steps) 
                        next_interval_state_idx = actual_state_sequence_zvs(k_interval + 1);
                        next_interval_switch_config = prms_zvs.valid_states(next_interval_state_idx, :);
                    else 
                        if ~isempty(actual_state_sequence_zvs) 
                            next_interval_state_idx = actual_state_sequence_zvs(1); 
                            next_interval_switch_config = prms_zvs.valid_states(next_interval_state_idx, :);
                        else
                            continue; 
                        end
                    end
                    
                    is_switch_off_in_current = (current_interval_switch_config(j_sw) == 0);
                    is_switch_on_in_next = (next_interval_switch_config(j_sw) == 1);
                    
                    is_leg_specific_deadtime = false;
                    complementary_idx_for_j_sw = NaN;
                    if isKey(complementary_map_local, j_sw)
                        complementary_idx_for_j_sw = complementary_map_local(j_sw);
                    end

                    if ~isnan(complementary_idx_for_j_sw) && complementary_idx_for_j_sw > 0 && complementary_idx_for_j_sw <= prms_zvs.num_switches
                        if current_interval_switch_config(j_sw) == 0 && current_interval_switch_config(complementary_idx_for_j_sw) == 0
                            is_leg_specific_deadtime = true;
                        end
                    else
                        if ismember(current_interval_state_idx, prms_zvs.topology_info.deadtime_states)
                            is_leg_specific_deadtime = true; 
                        end
                    end
                    
                    if is_leg_specific_deadtime && is_switch_off_in_current && is_switch_on_in_next
                        x_end_target_dt_interval = state_vectors_at_interval_ends{k_interval};
                        
                        if isempty(x_end_target_dt_interval) || target_cap_var_idx > length(x_end_target_dt_interval)
                             if ~prms_zvs.suppress_console_output
                                fprintf('    simulateZVSConditions: ERROR - x_end_target_dt_interval invalid or target_cap_var_idx out of bounds for Switch %d.\n', j_sw);
                             end
                             % Use adaptable penalty based on primary input voltage
                             penalty_V_err = 100.0; % Default if input voltage cannot be determined
                             if ~isempty(prms_zvs.u_vars) && ~isempty(prms_zvs.input_values) && isKey(prms_zvs.input_values, prms_zvs.u_vars{1})
                                 val_err = prms_zvs.input_values(prms_zvs.u_vars{1});
                                 if isscalar(val_err) && isnumeric(val_err) && abs(val_err) > 1e-9 % Ensure it's a usable number
                                     penalty_V_err = abs(val_err); 
                                 end
                             end
                             capacitor_voltages_per_switch(j_sw) = penalty_V_err; 
                             voltage_rates_per_switch(j_sw) = penalty_V_err / (0.05 * prms_zvs.T_period);
                        else
                            capacitor_voltages_per_switch(j_sw) = x_end_target_dt_interval(target_cap_var_idx);
                            
                            A_dt_interval = As{current_interval_state_idx};
                            B_dt_interval = Bs{current_interval_state_idx};
                            
                            dxdt_at_end_of_dt = A_dt_interval * x_end_target_dt_interval + B_dt_interval * u_vec_zvs;
                            voltage_rates_per_switch(j_sw) = dxdt_at_end_of_dt(target_cap_var_idx);
                        end
                        
                        found_zvs_interval_for_switch = true;
                        
                        if j_sw == 1 && ~prms_zvs.suppress_console_output 
                            switch_name_diag = prms_zvs.switches_ordered_list(j_sw).name;
                            fprintf('  DIAGNOSTIC (simulateZVSConditions for %s, j_sw=%d):\n', switch_name_diag, j_sw);
                            fprintf('    Dynamic Interval k_interval = %d (State %d -> %d)\n', k_interval, current_interval_state_idx, next_interval_state_idx);
                            fprintf('    Duration of this deadtime interval: %.4e s\n', actual_interval_durations_zvs(k_interval));
                            fprintf('    ZVS Vc(%s_Cs_parasitic) = %.6e V\n', switch_name_diag, capacitor_voltages_per_switch(j_sw));
                            fprintf('    ZVS dVcdt(%s_Cs_parasitic) = %.6e V/s\n', switch_name_diag, voltage_rates_per_switch(j_sw));
                        end
                        
                        if ~prms_zvs.suppress_console_output
                            fprintf('    simulateZVSConditions: ZVS condition for Switch %d (%s) evaluated at end of dynamic interval %d (State %d -> %d).\n', ...
                                    j_sw, prms_zvs.switches_ordered_list(j_sw).name, k_interval, current_interval_state_idx, next_interval_state_idx);
                        end
                        break; 
                    end
                end % End loop k_interval

                if ~found_zvs_interval_for_switch
                    if ~prms_zvs.suppress_console_output
                        fprintf('    simulateZVSConditions: NO ZVS INTERVAL for Switch %d (%s). Assigning ROBUST PENALTY.\n', j_sw, prms_zvs.switches_ordered_list(j_sw).name);
                    end
                    penalty_V_no_interval = 100.0; 
                     if ~isempty(prms_zvs.u_vars) && ~isempty(prms_zvs.input_values) && isKey(prms_zvs.input_values, prms_zvs.u_vars{1})
                        val_no_int = prms_zvs.input_values(prms_zvs.u_vars{1});
                        if isscalar(val_no_int) && isnumeric(val_no_int) && abs(val_no_int) > 1e-9
                            penalty_V_no_interval = abs(val_no_int) * 2; 
                        end
                    end
                    capacitor_voltages_per_switch(j_sw) = penalty_V_no_interval; 
                    voltage_rates_per_switch(j_sw) = penalty_V_no_interval / (0.05 * prms_zvs.T_period); 
                     if ~prms_zvs.suppress_console_output
                        fprintf('    simulateZVSConditions: PENALTY (No ZVS Interval) for Switch %d (%s): Vc=%.4e, dVcdt=%.4e\n', ...
                            j_sw, prms_zvs.switches_ordered_list(j_sw).name, ...
                            capacitor_voltages_per_switch(j_sw), voltage_rates_per_switch(j_sw));
                    end
                end
            end 
        end 
        
        if ~prms_zvs.suppress_console_output
            fprintf('    simulateZVSConditions: Finished. Voltages: [%s], Rates: [%s]\n', mat2str(capacitor_voltages_per_switch,4), mat2str(voltage_rates_per_switch,4));
        end
    end % end simulateZVSConditions
%-------------------------------------------------------------------------

    % Helper function to simulate one cycle of a dynamic sequence and get states at end of each interval
    function state_vectors_at_interval_ends = simulate_dynamic_sequence_step_ends(current_prms, x0_cycle_start, As_matrices, Bs_matrices, sequence_to_sim, durations_to_sim)
        % Simulates one cycle of the given dynamic sequence and durations.
        % Returns a cell array of state vectors at the end of each interval.
        
        num_dynamic_steps = length(sequence_to_sim);
        state_dim_local_sim = current_prms.state_dim;
        u_vec_sim = constructInputVector(current_prms);

        state_vectors_at_interval_ends = cell(1, num_dynamic_steps);
        current_state_vec_sim = x0_cycle_start; % Start with the provided initial condition

        for k_sim_step = 1:num_dynamic_steps
            state_idx_sim = sequence_to_sim(k_sim_step);
            dt_interval_sim = durations_to_sim(k_sim_step);

            if state_idx_sim < 1 || state_idx_sim > length(As_matrices) || isempty(As_matrices{state_idx_sim}) || ...
               state_idx_sim > length(Bs_matrices) || isempty(Bs_matrices{state_idx_sim})
                error('simulate_dynamic_sequence:MissingStateModel', 'State model missing for state %d in dynamic sequence simulation.', state_idx_sim);
            end

            A_sim = As_matrices{state_idx_sim};
            B_sim = Bs_matrices{state_idx_sim};

            if dt_interval_sim > 1e-20 % Propagate only if duration is significant
                lambda_k_sim = calculateLambda(A_sim, dt_interval_sim, current_prms);
                Phi_k_sim = eye(state_dim_local_sim) + A_sim * lambda_k_sim;
                
                if size(B_sim, 2) == current_prms.num_inputs
                    Gamma_k_sim = lambda_k_sim * B_sim * u_vec_sim;
                    current_state_vec_sim = Phi_k_sim * current_state_vec_sim + Gamma_k_sim;
                else
                    error('simulate_dynamic_sequence:DimensionMismatchB', 'Dimension mismatch B matrix for state %d in dynamic sequence.', state_idx_sim);
                end
            end
            % Store the state vector at the end of this interval
            state_vectors_at_interval_ends{k_sim_step} = current_state_vec_sim;
        end
    end % end simulate_dynamic_sequence_step_ends


%------------------------------------------------------------------
function [optimal_dt_out, dt_history_out, voltage_history_out] = findOptimalDeadtimes(prms, As, Bs)
        % findOptimalDeadtimes: Optimizes per-switch applied deadtimes for ZVS.
        % Uses event-driven timing method for sequence generation.
        % Applies stricter individual deadtime bounding and configurable max step change.

        % --- Extract parameters ---
        max_iterations = prms.max_iterations;
        voltage_tolerance = prms.voltage_tolerance;
        damping_factor = prms.damping_factor;
        min_dt_individual_switch = prms.min_dt; 
        
        T_period_local = prms.T_period; 
        
        % Absolute maximum for an individual deadtime
        if isfield(prms, 'max_individual_dt_fraction_of_T')
            max_individual_dt_val = T_period_local * prms.max_individual_dt_fraction_of_T;
        else
            % This warning should ideally be caught by validateAndMergeConfig if default is used
            warning('FindOptimalDeadtimes:MissingMaxIndivDTFraction', 'prms.max_individual_dt_fraction_of_T not found. Defaulting to 0.25 * T_period for max individual deadtime.');
            max_individual_dt_val = T_period_local * 0.25; % Default if somehow missed by validation
        end
        
        % Maximum change allowed for a deadtime in one iteration, calculated from the fraction in prms
        % The prms.max_dt_step_change_fraction is validated in validateAndMergeConfig
        max_dt_step_change = T_period_local * prms.max_dt_step_change_fraction; 

        num_switches_total = prms.num_switches;

        % --- Initialize Optimization Variables ---
        if length(prms.dt) ~= num_switches_total
            error('FindOptimalDeadtimes:InitialDtLengthMismatch', ...
                  'Length of prms.dt (initial per-switch deadtimes) (%d) must match num_switches (%d).', length(prms.dt), num_switches_total);
        end
        dt_current_per_switch = prms.dt(:)'; 
        
        dt_history_out = zeros(max_iterations, num_switches_total);
        voltage_history_out = NaN(max_iterations, num_switches_total); 

        if ~prms.suppress_console_output
            fprintf('  Starting Newton-Raphson for ZVS deadtime optimization (%d switches)...\n', num_switches_total);
            fprintf('    Max individual deadtime: %.3f ns, Max step change: %.3f ns (%.1f%% of T_period)\n', ...
                max_individual_dt_val*1e9, max_dt_step_change*1e9, prms.max_dt_step_change_fraction*100);
            dt_header_fmt_str = repmat('%-9.3f ', 1, num_switches_total);
            vc_header_fmt_str = repmat('%-11.4g ', 1, num_switches_total);
            header_part1 = sprintf('  %4s | %-*s |', 'Iter', length(sprintf(dt_header_fmt_str, ones(1,num_switches_total)*1e9)), 'Applied Deadtimes (ns) per Switch');
            header_part2 = sprintf(' %-*s | %-15s\n', length(sprintf(vc_header_fmt_str, ones(1,num_switches_total))), 'Target Voltages (V) per Switch', 'Max |dV/dt|');
            fprintf('%s%s', header_part1, header_part2);
            underline_part1 = sprintf('  %4s | %s |', '----', repmat('-',1,length(sprintf(dt_header_fmt_str, ones(1,num_switches_total)*1e9))));
            underline_part2 = sprintf(' %s | %s\n', repmat('-',1,length(sprintf(vc_header_fmt_str, ones(1,num_switches_total)))), repmat('-',1,15) );
            fprintf('%s%s', underline_part1, underline_part2);
        end

        iter_params = prms; 
        converged = false;
        iter_count = 0; 

        is_switch_targeted_for_zvs = false(1, num_switches_total);
        zvs_cap_state_var_indices = NaN(1, num_switches_total); 

        if isfield(prms, 'topology_info') && isfield(prms.topology_info, 'zvs_target_map_switch_to_cap_idx') && ...
           isa(prms.topology_info.zvs_target_map_switch_to_cap_idx, 'containers.Map')
            map_sw_idx_to_cap_var_idx = prms.topology_info.zvs_target_map_switch_to_cap_idx;
            for target_sw_name_cell = prms.zvs_switches_to_target 
                target_sw_name = target_sw_name_cell{1};
                found_sw_ordered_idx = 0;
                for k_sw_ordered = 1:num_switches_total
                    if strcmp(prms.switches_ordered_list(k_sw_ordered).name, target_sw_name)
                        found_sw_ordered_idx = k_sw_ordered;
                        break;
                    end
                end
                if found_sw_ordered_idx > 0
                    if isKey(map_sw_idx_to_cap_var_idx, found_sw_ordered_idx)
                        cap_var_idx = map_sw_idx_to_cap_var_idx(found_sw_ordered_idx);
                        if ~isnan(cap_var_idx) && cap_var_idx > 0 && cap_var_idx <= prms.state_dim
                            is_switch_targeted_for_zvs(found_sw_ordered_idx) = true;
                            zvs_cap_state_var_indices(found_sw_ordered_idx) = cap_var_idx;
                            if ~prms.suppress_console_output; fprintf('    Targeting ZVS for switch %s (idx %d), ZVS Cap Var Idx: %d.\n', target_sw_name, found_sw_ordered_idx, cap_var_idx); end
                        else
                             if ~prms.suppress_console_output; warning('FindOptimalDeadtimes:InvalidCapVarIdxInMap', 'Switch %s (idx %d) has invalid ZVS cap var idx (%g) in map.', target_sw_name, found_sw_ordered_idx, cap_var_idx); end
                        end
                    else
                         if ~prms.suppress_console_output; warning('FindOptimalDeadtimes:SwitchNotInMap', 'Switch %s (idx %d) targeted by user for ZVS not found in zvs_target_map_switch_to_cap_idx.',target_sw_name, found_sw_ordered_idx); end
                    end
                else
                    if ~prms.suppress_console_output; warning('FindOptimalDeadtimes:TargetSwitchNameNotFound', 'User-targeted ZVS switch name "%s" not found in prms.switches_ordered_list.', target_sw_name); end
                end
            end
        else
            if ~prms.suppress_console_output; warning('FindOptimalDeadtimes:MissingZVSTargetMapStruct', 'prms.topology_info.zvs_target_map_switch_to_cap_idx is missing. Cannot determine ZVS targets.'); end
            optimal_dt_out = dt_current_per_switch; dt_history_out = []; voltage_history_out = [];
            if prms.run_optimization && ~prms.suppress_console_output; fprintf('  ZVS Optimization disabled due to missing target map.\n'); end
            return;
        end
        
        num_actually_targeted_switches = sum(is_switch_targeted_for_zvs);
        if num_actually_targeted_switches == 0 && prms.run_optimization
            if ~prms.suppress_console_output; fprintf('  No switches were successfully mapped for ZVS targeting. Skipping optimization loop.\n'); end
            optimal_dt_out = dt_current_per_switch; dt_history_out = []; voltage_history_out = [];
            return;
        end

        for iter_count = 1:max_iterations
            dt_history_out(iter_count, :) = dt_current_per_switch;
            
            iter_params.dt = dt_current_per_switch; 
            iter_params.deadtimes_applied = dt_current_per_switch; 

            if prms.state_dim > 0
                 x0_iter = enhancedSVA(iter_params, As, Bs); 
            else
                 x0_iter = zeros(0,1); 
            end
            
            iter_params.is_switch_targeted_for_zvs = is_switch_targeted_for_zvs;
            iter_params.zvs_cap_state_var_indices = zvs_cap_state_var_indices;
            
            [target_voltages_iter, target_rates_iter] = simulateZVSConditions(iter_params, x0_iter, As, Bs); 

            voltage_history_out(iter_count, :) = target_voltages_iter;

            if any(isnan(target_voltages_iter(is_switch_targeted_for_zvs))) || any(isnan(target_rates_iter(is_switch_targeted_for_zvs)))
                 nan_indices = find( (isnan(target_voltages_iter) | isnan(target_rates_iter)) & is_switch_targeted_for_zvs );
                 if ~prms.suppress_console_output 
                    warning('AnalyzeSwitchingConverter:ZVSsimNaN', 'ZVS sim returned NaN for targeted switches (indices: %s) at iter %d. Stopping.', mat2str(nan_indices), iter_count);
                 end
                 iter_count = iter_count - 1; 
                 break; 
            end

            if ~prms.suppress_console_output
                dt_str_disp = sprintf(dt_header_fmt_str, dt_current_per_switch*1e9);
                vc_str_disp = sprintf(vc_header_fmt_str, target_voltages_iter);
                valid_rates_for_disp = target_rates_iter(is_switch_targeted_for_zvs & ~isnan(target_rates_iter));
                max_abs_rate_disp = 0;
                if ~isempty(valid_rates_for_disp); max_abs_rate_disp = max(abs(valid_rates_for_disp*1e-9)); end 
                fprintf('  %4d | %s | %s | %15.3g\n', iter_count, dt_str_disp, vc_str_disp, max_abs_rate_disp);
            end

            if num_actually_targeted_switches > 0
                converged_flags = abs(target_voltages_iter(is_switch_targeted_for_zvs)) < voltage_tolerance;
                if all(converged_flags)
                    if ~prms.suppress_console_output; fprintf('  Convergence achieved for all %d targeted switch(es) after %d iterations.\n', num_actually_targeted_switches, iter_count); end
                    converged = true;
                    break;
                end
            else 
                converged = true; 
                break; 
            end

            dt_new_per_switch = dt_current_per_switch; 
            for j_sw = 1:num_switches_total
                if is_switch_targeted_for_zvs(j_sw)
                    Vc_j = target_voltages_iter(j_sw);
                    dVcdt_j = target_rates_iter(j_sw);
                    dt_update_step = 0; 

                    if abs(dVcdt_j) > prms.significant_rate 
                        dt_update_step = -damping_factor * (Vc_j / dVcdt_j);
                    else 
                        if iter_count > 1 && ~prms.suppress_console_output && abs(Vc_j) > voltage_tolerance 
                             warning('FindOptimalDeadtimes:NearZeroRateHeuristic', 'Near-zero rate (%.2e V/s) for sw %d (%s) at iter %d with Vc=%.2f V. Heuristic used.', dVcdt_j, j_sw, prms.switches_ordered_list(j_sw).name, iter_count, Vc_j);
                        end
                        current_dt_val = dt_current_per_switch(j_sw);
                        if Vc_j > voltage_tolerance; dt_update_step = current_dt_val * (prms.adjust_increase - 1); 
                        elseif Vc_j < -voltage_tolerance; dt_update_step = current_dt_val * (prms.adjust_decrease - 1); 
                        end
                    end
                    
                    % Apply the maximum step change limit
                    if abs(dt_update_step) > max_dt_step_change
                        if ~prms.suppress_console_output && max_dt_step_change < T_period_local % Only print if limit is meaningfully restrictive
                            fprintf('    Limiting dt change for switch %d from %.3e to +/- %.3e ns.\n', j_sw, dt_update_step*1e9, max_dt_step_change*1e9);
                        end
                        dt_update_step = sign(dt_update_step) * max_dt_step_change;
                    end
                    dt_new_per_switch(j_sw) = dt_current_per_switch(j_sw) + dt_update_step;
                end 
            end 

            dt_new_per_switch = max(min_dt_individual_switch, dt_new_per_switch); 
            dt_new_per_switch = min(max_individual_dt_val, dt_new_per_switch);
            
            dt_current_per_switch = dt_new_per_switch;
        end 

        if ~converged 
            if iter_count == max_iterations && ~prms.suppress_console_output
                fprintf('  Warning: Max iterations (%d) reached for ZVS optimization without full convergence.\n', max_iterations);
            elseif iter_count < max_iterations && iter_count >= 0 
                 if ~prms.suppress_console_output && iter_count > 0 
                    fprintf('  Warning: ZVS Optimization stopped prematurely at iteration %d.\n', iter_count); 
                 elseif ~prms.suppress_console_output && iter_count == 0
                    fprintf('  Warning: ZVS Optimization stopped before completing the first iteration.\n');
                 end
            end
            
            if iter_count > 0 
                voltages_of_targeted_switches_hist = voltage_history_out(1:iter_count, is_switch_targeted_for_zvs);
                
                if ~isempty(voltages_of_targeted_switches_hist) && any(~isnan(voltages_of_targeted_switches_hist(:)))
                    sum_abs_errors_hist = sum(abs(voltages_of_targeted_switches_hist), 2, 'omitnan');
                    
                    valid_error_sum_indices = find(~isnan(sum_abs_errors_hist));
                    
                    if ~isempty(valid_error_sum_indices)
                        [~, min_idx_within_valid] = min(sum_abs_errors_hist(valid_error_sum_indices));
                        best_idx = valid_error_sum_indices(min_idx_within_valid);
                        
                        dt_current_per_switch = dt_history_out(best_idx, :);
                        if ~prms.suppress_console_output
                            fprintf('  Using best deadtimes from iteration %d: [%s] ns (Error sum for targeted: %.4g)\n', ...
                                    best_idx, sprintf('%.3f ', dt_current_per_switch*1e9), sum_abs_errors_hist(best_idx));
                        end
                    else
                         if ~prms.suppress_console_output
                            fprintf('  Could not determine best iteration (all error sums NaN for targeted switches). Using deadtimes from last valid iteration %d.\n', iter_count);
                         end
                         dt_current_per_switch = dt_history_out(iter_count, :); 
                    end
                else
                     if ~prms.suppress_console_output
                        fprintf('  No valid voltage history for any targeted switches. Using deadtimes from last valid iteration %d.\n', iter_count);
                     end
                     dt_current_per_switch = dt_history_out(iter_count, :);
                end
            else 
                if ~prms.suppress_console_output
                    fprintf('  Optimization did not complete any iterations with history. Returning initial deadtimes.\n');
                end
                dt_current_per_switch = prms.dt(:)'; 
            end
        end

        optimal_dt_out = dt_current_per_switch;
        
        if iter_count > 0
            dt_history_out = dt_history_out(1:iter_count, :);
            voltage_history_out = voltage_history_out(1:iter_count, :);
        else 
            dt_history_out = []; 
            voltage_history_out = [];
        end

    end % end findOptimalDeadtimes



   %% --- Phase 5 Helpers ---
 function x0 = enhancedSVA(prms, As, Bs)
        % enhancedSVA: Calculates steady-state initial conditions using ESVA.
        % Uses event-driven sequence and durations from generateSequenceFromTimings.
        %
        % Inputs:
        %   prms - Struct containing all necessary parameters:
        %     .ts_switches          - Cell array: Ideal ON/OFF timings for each switch.
        %     .deadtimes_applied    - Vector: Current deadtimes for ON transitions, per switch.
        %                             (This is expected to be prms.dt from the calling context).
        %     .switches_ordered_list - Struct array: Info about switches.
        %     .T_period             - Scalar: Switching period.
        %     .valid_states         - Matrix: Valid switch state vectors.
        %     .state_dim            - Scalar: Dimension of the state vector.
        %     .num_valid_states     - Scalar: Total number of unique valid states.
        %     .num_inputs           - Scalar: Number of input sources.
        %     .abstol               - Scalar: Absolute tolerance for calculateLambda.
        %     (and other fields needed by generateSequenceFromTimings and calculateLambda)
        %
        %   As     - Cell array: A matrices for each valid state.
        %   Bs     - Cell array: B matrices for each valid state.
        %
        % Output:
        %   x0     - Column vector: Steady-state initial condition.

        % Handle case with no state variables
        if prms.state_dim == 0
            x0 = zeros(0,1); % Return empty column vector
            if ~prms.suppress_console_output
                fprintf('    enhancedSVA: No state variables, returning empty x0.\n');
            end
            return;
        end

        % --- Step 1: Generate State Sequence and Durations using Event-Driven Method ---
        % The 'prms' struct should contain all necessary fields for generateSequenceFromTimings.
        % prms.deadtimes_applied should hold the current per-switch deadtimes.
        
        prms_for_seq_gen = prms; % Pass a copy of prms
        % Ensure deadtimes_applied is correctly set if prms.dt is the source
        if isfield(prms, 'dt') && ~isfield(prms, 'deadtimes_applied')
            prms_for_seq_gen.deadtimes_applied = prms.dt;
        elseif ~isfield(prms, 'deadtimes_applied') && isfield(prms, 'initial_deadtimes')
            % Fallback if .dt isn't set but .initial_deadtimes (per-switch) is
            prms_for_seq_gen.deadtimes_applied = prms.initial_deadtimes;
        elseif ~isfield(prms_for_seq_gen, 'deadtimes_applied')
             error('enhancedSVA:MissingDeadtimesApplied', 'prms.deadtimes_applied is missing for sequence generation.');
        end


        if ~prms.suppress_console_output
            fprintf('    enhancedSVA: Calling generateSequenceFromTimings (event-driven method)...\n');
        end
        [actual_state_sequence, actual_interval_durations] = generateSequenceFromTimings(prms_for_seq_gen); % Nested function call

        if isempty(actual_state_sequence)
            error('enhancedSVA:EmptySequence', 'generateSequenceFromTimings (event-driven method) returned an empty sequence. Cannot proceed.');
        end
        if ~prms.suppress_console_output
            fprintf('    enhancedSVA: Generated sequence (event-driven) length: %d\n', length(actual_state_sequence));
        end

        num_seq_steps = length(actual_state_sequence);
        state_dim_local = prms.state_dim;

        % --- Step 2: Pre-calculate Phi_k and Gamma_k for each step in the generated sequence ---
        Phis_seq = cell(1, num_seq_steps); 
        Gammas_seq = cell(1, num_seq_steps); 

        u_vec = constructInputVector(prms); % Construct the input vector u ONCE

        for k_esva = 1:num_seq_steps
            state_idx = actual_state_sequence(k_esva);    % Get the actual state index from generated sequence
            dt_interval = actual_interval_durations(k_esva); % Get the duration from generated durations

            if state_idx < 1 || state_idx > prms.num_valid_states || ...
               state_idx > length(As) || isempty(As{state_idx}) || ...
               state_idx > length(Bs) || isempty(Bs{state_idx})
                 error('ESVA:MissingStateModel', 'State model missing for state %d (step %d in generated sequence).', state_idx, k_esva);
            end
            A = As{state_idx}; 
            B = Bs{state_idx};

            if dt_interval > 1e-20 % Use a small threshold to avoid issues with zero duration
                lambda_k = calculateLambda(A, dt_interval, prms); 
                Phis_seq{k_esva} = eye(state_dim_local) + A * lambda_k; 
                if size(B, 2) == prms.num_inputs 
                    Gammas_seq{k_esva} = lambda_k * B * u_vec;
                else
                     error('ESVA:DimensionMismatchB', 'Dimension mismatch between B matrix (%d cols) and number of inputs (%d) for state %d.', size(B,2), prms.num_inputs, state_idx);
                end
            else 
                Phis_seq{k_esva} = eye(state_dim_local); 
                Gammas_seq{k_esva} = zeros(state_dim_local, 1); 
            end
        end 

        % --- Step 3: Calculate composite state transition matrix over one period ---
        Phi_composite = eye(state_dim_local);
        for k_phi_comp = 1:num_seq_steps 
             Phi_k_current = Phis_seq{k_phi_comp}; 
             Phi_composite = Phi_k_current * Phi_composite; 
        end

        % --- Step 4: Calculate composite input effect term over one period ---
        Gamma_composite = zeros(state_dim_local, 1);
        Phi_partial_from_k_plus_1_to_end = eye(state_dim_local); 
        for k_gamma_comp = num_seq_steps:-1:1 
             Gamma_k = Gammas_seq{k_gamma_comp}; 
             Gamma_composite = Gamma_composite + Phi_partial_from_k_plus_1_to_end * Gamma_k;
             Phi_k_stored = Phis_seq{k_gamma_comp}; 
             Phi_partial_from_k_plus_1_to_end = Phi_partial_from_k_plus_1_to_end * Phi_k_stored;
        end

        % --- Step 5: Solve for steady-state initial condition ---
        I_mat = eye(state_dim_local); 
        matrix_to_invert = (I_mat - Phi_composite);
        
        matrix_cond = 0; % Default if N_mna is 0 or matrix_to_invert is empty
        if ~isempty(matrix_to_invert)
            matrix_cond = rcond(matrix_to_invert);
        end

        if ~prms.suppress_console_output
            fprintf('    enhancedSVA: Condition number of (I - Phi_composite): %.3e\n', matrix_cond);
            if matrix_cond < prms.degeneracy_threshold 
                 warning('AnalyzeSwitchingConverter:ESVASingularMatrix', ...
                      'Matrix (I - Phi_composite) is singular or ill-conditioned (rcond=%.2e < threshold %.2e). PSS solution might be inaccurate. Using pseudo-inverse.', matrix_cond, prms.degeneracy_threshold);
            end
        end
        
        if isempty(matrix_to_invert) % Should only happen if state_dim_local is 0, handled at start
            x0 = zeros(0,1);
        elseif rcond(matrix_to_invert) >= prms.degeneracy_threshold 
            x0 = matrix_to_invert \ Gamma_composite;
        else
            if ~prms.suppress_console_output 
                 warning('AnalyzeSwitchingConverter:ESVASingularMatrixFallback', ...
                      'Matrix (I - Phi_composite) rcond is very low (%.2e). Using pseudo-inverse (pinv).', rcond(matrix_to_invert));
            end
            x0 = pinv(matrix_to_invert) * Gamma_composite;
        end
        if ~prms.suppress_console_output
            fprintf('    enhancedSVA: x0 calculated.\n');
        end

    end % end enhancedSVA


    %---------------------------------------------
 function lambda = calculateLambda(A, dt, prms) % <<< MODIFIED: Added prms argument
        % Calculates the lambda term: integral(expm(A*t), t=0..dt)
        % using Taylor series with adaptive order selection and scaling/squaring.
        % Based on Sadri (2022) [cite: 33-35] and Ellouz_ESVA_Analysis_NetlistImport_Update9_4.txt [cite: 637-688].

        % Handle edge cases: zero duration or empty matrix
        if dt <= 1e-20 % Treat very small dt as zero
            lambda = zeros(size(A));
            return;
        end
        if isempty(A) % Handle case where A might be empty (e.g., state_dim=0)
             lambda = [];
             return;
        end

        % --- Scaling (similar to expm, prevents issues with large A*dt norm) ---
        normA = norm(A, 'inf'); % Estimate matrix norm
        scaling = 0;
        % Heuristic: scale if norm(A*dt) is large enough (e.g., > 0.5 or 1)
        norm_A_dt = normA * dt;
        if norm_A_dt > 0.5
            % Determine scaling factor s such that norm(A*dt / 2^s) is small
            scaling = max(0, ceil(log2(norm_A_dt)));
        end
        dt_scaled = dt / (2^scaling); % Apply scaling

        % --- Taylor Series Calculation (for scaled dt) ---
        % lambda_approx = sum_{n=1 to k} A^(n-1) * dt_scaled^n / n!
        % Calculate iteratively: term(n) = term(n-1) * [ A * dt_scaled / n ]
        I = eye(size(A));
        A_dt_scaled_term = A * dt_scaled; % Pre-calculate A * dt_scaled

        term_k = I * dt_scaled; % First term (n=1): I * dt_scaled^1 / 1!
        lambda_approx = term_k; % Initialize sum with the first term
        k = 1; % Current term order (n)
        max_terms = 100; % Limit iterations to prevent infinite loops

        % <<< MODIFIED: Use prms argument to get tolerance >>>
        if ~isfield(prms, 'abstol')
             error('InternalError:MissingAbstol', 'abstol field missing from params struct passed to calculateLambda.');
        end
        abstol_lambda = prms.abstol;
        % <<< END MODIFICATION >>>

        while k < max_terms
            % Calculate next term (n=k+1) using previous term (n=k)
            % term(k+1) = term(k) * [ A * dt_scaled / (k+1) ]
            term_k_plus_1 = term_k * (A_dt_scaled_term / (k + 1));
            lambda_new = lambda_approx + term_k_plus_1; % Update sum

            % Check convergence based on the norm of the added term relative to current sum
            norm_term = norm(term_k_plus_1, 'inf');
            norm_lambda = norm(lambda_new, 'inf');
            % Use relative tolerance, but handle case where norm_lambda is near zero
            rel_change = ifelse(norm_lambda > eps, norm_term / norm_lambda, norm_term);

            lambda_approx = lambda_new; % Store updated approximation
            term_k = term_k_plus_1; % Store current term for next iteration's calculation

            % Break if change is below tolerance (and we've done at least a few terms)
            if rel_change < abstol_lambda && k > 2
                break;
            end
            k = k + 1; % Increment term order
        end % end while Taylor series

        if k == max_terms
             warning('AnalyzeSwitchingConverter:LambdaMaxTerms', 'Lambda calculation reached max terms (%d) for dt=%.2e. Result might be inaccurate.', max_terms, dt);
        end

        % --- Undo scaling using the squaring phase relation ---
        % lambda(2t) = (expm(A*t) + I) * lambda(t)
        % where expm(A*t) = A*lambda(t) + I (approximately from Taylor series relation)
        % => lambda(2t) = (A*lambda(t) + 2*I) * lambda(t) = 2*lambda(t) + A*lambda(t)*lambda(t)
        lambda = lambda_approx; % Start with lambda for dt_scaled
        if scaling > 0
            for s_iter = 1:scaling
                 % Apply the doubling formula s times
                 lambda = 2*lambda + A*lambda*lambda;
            end
        end

    end % end calculateLambda

    %------------------------------------------------
    function intervals = calculateIntervals(prms)
        % Calculates the actual time durations for each unique state index
        % for the *current* simulation step (using current prms.dt).
        % It distributes the available non-deadtime duration according to
        % the proportions specified in prms.state_duration_proportions.
        % Sets intervals for deadtime states to ZERO.

        % --- Get necessary parameters ---
        T = 1/prms.fs; % Switching period
        num_states_local = prms.num_valid_states; % Total number of unique valid states
        state_sequence = prms.state_sequence; % The sequence of state indices
        proportions_map = prms.state_duration_proportions; % Map {'state_idx_str'} -> proportion
        current_deadtimes = prms.dt(:)'; % Current deadtime vector (ensure row)

        deadtime_states_indices = []; % Indices of states considered deadtime
        deadtime_end_seq_indices = []; % Indices in sequence where deadtimes end
        if isfield(prms, 'topology_info')
            if isfield(prms.topology_info, 'deadtime_states'); deadtime_states_indices = prms.topology_info.deadtime_states; end
            if isfield(prms.topology_info, 'deadtime_end_seq_indices'); deadtime_end_seq_indices = prms.topology_info.deadtime_end_seq_indices; end
        end
        num_deadtime_occurrences = length(deadtime_end_seq_indices);

        % --- Calculate total current deadtime based on sequence ---
        total_current_dt = 0;
        if num_deadtime_occurrences > 0
            if length(current_deadtimes) == num_deadtime_occurrences
                total_current_dt = sum(current_deadtimes);
            else
                % This mismatch should ideally be caught earlier
                warning('CalculateIntervals:DtLengthMismatch', 'Length of current prms.dt (%d) mismatches number of deadtime occurrences (%d). Assuming zero total deadtime.', length(current_deadtimes), num_deadtime_occurrences);
                total_current_dt = 0;
            end
        end
        total_current_dt = max(0, total_current_dt); % Ensure non-negative

        % --- Calculate time available for non-deadtime states ---
        T_available_for_non_deadtime = T - total_current_dt;
        if T_available_for_non_deadtime < -T*1e-9 % Allow small negative tolerance
            warning('CalculateIntervals:NegativeAvailableTime', ...
                  'Current total deadtime (%.4gus) exceeds period T (%.4gus). Setting non-deadtime state durations to zero.', ...
                  total_current_dt*1e6, T*1e6);
            T_available_for_non_deadtime = 0;
        elseif T_available_for_non_deadtime < 0
            T_available_for_non_deadtime = 0; % Clamp to zero
        end

        % --- Initialize the output array ---
        intervals = zeros(1, num_states_local);

        % --- Distribute available time to non-deadtime states based on proportions ---
        non_deadtime_indices_in_map = cellfun(@str2double, keys(proportions_map));

        for i = non_deadtime_indices_in_map
             % Check if state index is valid and not a deadtime state
             if i > 0 && i <= num_states_local && ~ismember(i, deadtime_states_indices)
                 state_idx_str = num2str(i);
                 if isKey(proportions_map, state_idx_str)
                     proportion = proportions_map(state_idx_str);
                     intervals(i) = proportion * T_available_for_non_deadtime;
                     % Ensure non-negative due to potential floating point issues
                     intervals(i) = max(0, intervals(i));
                 else
                      % This should be caught by validation if the state is in the sequence
                      if ismember(i, state_sequence)
                          warning('CalculateIntervals:MissingProportion', 'Proportion for non-deadtime state index %d not found in map, but state is in sequence. Assigning zero duration.', i);
                      end
                      intervals(i) = 0;
                 end
             elseif ismember(i, deadtime_states_indices)
                  % This case should not happen if map keys are validated correctly
                  warning('CalculateIntervals:DeadtimeInProportionMap', 'Deadtime state index %d found as key in state_duration_proportions map. Ignoring.', i);
             end
        end

        % Deadtime states inherently have intervals(i) = 0 in this array

        % --- Final Sum Check (Optional Debugging/Sanity Check) ---
        % This check verifies the logic of *this* function.
        % final_sum_check = 0;
        % dt_check_idx = 1;
        % for k_seq_chk = 1:length(state_sequence)
        %     idx_chk = state_sequence(k_seq_chk);
        %     if idx_chk > 0 && idx_chk <= num_states_local
        %         if ismember(idx_chk, deadtime_states_indices)
        %             if dt_check_idx <= length(current_deadtimes); final_sum_check = final_sum_check + current_deadtimes(dt_check_idx); dt_check_idx = dt_check_idx + 1; end
        %         else
        %             final_sum_check = final_sum_check + intervals(idx_chk);
        %         end
        %     end
        % end
        % if abs(final_sum_check - T) > T*1e-9
        %      warning('CalculateIntervals:InternalSumCheckFail', 'Internal sum check failed: %.6f us vs T=%.6f us.', final_sum_check*1e6, T*1e6);
        % end
        % --- End Final Sum Check ---

    end % end calculateIntervals



%------------------------------------------------
    function [t_out, x_out, y_out] = generateOutputWaveforms(prms, x0_ss_in, As, Bs, Cs, Ds)
        % generateOutputWaveforms: Simulates steady-state periods and calculates outputs.
        % MODIFIED:
        % - Uses generateSequenceFromTimings with final per-switch deadtimes (prms.dt)
        %   to get the actual sequence and durations for one period.
        % - Simulates point-by-point using this dynamically generated sequence, repeated.
        % - The old calculateIntervals and fixed state_sequence logic for waveform generation is removed.

        T_period = prms.T_period; % Period from prms (already 1/fs)
        num_periods_to_sim = 3;
        total_sim_time = num_periods_to_sim * T_period;

        points_per_period_target = 2000; 
        num_points = round(num_periods_to_sim * points_per_period_target);
        % Ensure at least two points if total_sim_time is very small or zero
        if total_sim_time > 0 && num_points < 2 
            num_points = 2;
        elseif total_sim_time == 0 && num_points < 1 % Handle zero simulation time
            num_points = 1;
        end
        
        if num_points <= 1 && total_sim_time > 0
             warning('generateOutputWaveforms:LowNumPoints', 'Number of simulation points is very low (%d) for total_sim_time %.4e. Waveform resolution might be poor.', num_points, total_sim_time);
             if num_points < 2; num_points = 2; end % Force at least 2 for linspace if time > 0
        end

        if total_sim_time > 0
            t_out = linspace(0, total_sim_time, num_points);
        else
            t_out = 0; % Single point at t=0 if total_sim_time is zero
            num_points = 1;
        end
        
        dt_sim_step = 0;
        if num_points > 1
            dt_sim_step = t_out(2) - t_out(1);
        end

        state_dim_local = prms.state_dim;
        num_valid_states_local = prms.num_valid_states;

        num_outputs = 0;
        if iscell(prms.outputs)
            if isscalar(prms.outputs) && strcmp(prms.outputs{1}, 'all_states')
                num_outputs = prms.state_dim;
            else
                num_outputs = length(prms.outputs);
            end
        end

        x_out = zeros(state_dim_local, num_points);
        y_out = zeros(num_outputs, num_points);

        % --- Handle Static Case (No Switches or only 1 valid state) ---
        if prms.num_switches == 0 || num_valid_states_local == 1
            if ~prms.suppress_console_output
                fprintf('    generateOutputWaveforms: Generating waveform for static or single-state circuit...\n');
            end
            if state_dim_local > 0
                A = As{1}; B = Bs{1}; C = []; D = [];
                if num_outputs > 0; C = Cs{1}; D = Ds{1}; end
                
                u_vec = constructInputVector(prms);
                x_out(:, 1) = x0_ss_in;
                if num_outputs > 0 && ~isempty(C) && ~isempty(D)
                    if size(C,2) == state_dim_local && size(D,2) == prms.num_inputs
                        y_out(:, 1) = C * x0_ss_in + D * u_vec;
                    else
                         warning('generateOutputWaveforms:StaticCDMismatchT0', 'Static C/D matrix dimensions mismatch at t=0.');
                    end
                end

                current_state_vec_wf = x0_ss_in;
                for k_t = 2:num_points
                    if dt_sim_step > 0
                        lambda_step = calculateLambda(A, dt_sim_step, prms);
                        Phi_step = eye(state_dim_local) + A * lambda_step;
                        if size(B, 2) == prms.num_inputs
                            Gamma_step = lambda_step * B * u_vec;
                            current_state_vec_wf = Phi_step * current_state_vec_wf + Gamma_step;
                        else
                             error('generateOutputWaveforms:StaticBMismatchLoop', 'Static B matrix dimension mismatch in loop.');
                        end
                    end % else, state remains x0_ss_in if dt_sim_step is 0
                    x_out(:, k_t) = current_state_vec_wf;
                    if num_outputs > 0 && ~isempty(C) && ~isempty(D)
                        if size(C,2) == state_dim_local && size(D,2) == prms.num_inputs
                             y_out(:, k_t) = C * current_state_vec_wf + D * u_vec;
                        % else warning already issued
                        end
                    end
                end
            else % No switches AND no states
                 if ~prms.suppress_console_output; fprintf('    generateOutputWaveforms: Skipping waveform generation: No switches and no state variables.\n'); end
                 if num_outputs > 0 && ~isempty(Ds) && ~isempty(Ds{1})
                     D_t0 = Ds{1}; u_vec_static = constructInputVector(prms);
                     if size(D_t0, 2) == prms.num_inputs; y_static = D_t0 * u_vec_static;
                         if size(y_static,1) == num_outputs; y_out(:, :) = repmat(y_static, 1, num_points); end
                     end
                 end
            end
            return; 
        end

        % --- Dynamic Sequence Generation for Waveforms (Switching Case) ---
        if ~prms.suppress_console_output
            fprintf('    generateOutputWaveforms: Generating dynamic sequence for waveforms using final deadtimes.\n');
        end
        
        prms_for_waveform_seq_gen = prms;
        % prms.dt contains the final per-switch deadtimes from params_final.dt
        prms_for_waveform_seq_gen.deadtimes_applied = prms.dt; 

        [actual_state_sequence_wf, actual_interval_durations_wf] = generateSequenceFromTimings(prms_for_waveform_seq_gen);

        if isempty(actual_state_sequence_wf)
            warning('generateOutputWaveforms:EmptyDynamicSequenceWF', 'generateSequenceFromTimings returned empty sequence for waveforms. Waveforms might be incorrect.');
            % Fallback: use first valid state for the whole period if sequence generation fails
            actual_state_sequence_wf = 1; 
            actual_interval_durations_wf = T_period;
            if num_valid_states_local == 0 % Should not happen if this point is reached
                x_out = NaN(state_dim_local, num_points); y_out = NaN(num_outputs, num_points); return;
            end
        end
        
        num_steps_in_dynamic_sequence = length(actual_state_sequence_wf);
        if ~prms.suppress_console_output
            fprintf('    Dynamic sequence for waveforms (length %d): %s\n', num_steps_in_dynamic_sequence, mat2str(actual_state_sequence_wf));
            fprintf('    Dynamic durations (us): %s\n', sprintf('%.4f ', actual_interval_durations_wf*1e6));
        end

        % --- Simulation Setup (Switching Case) ---
        x_out(:, 1) = x0_ss_in;
        current_state_vec_wf = x0_ss_in;
        u_vec = constructInputVector(prms);

        % Calculate output at t=0 using the first state of the dynamic sequence
        if num_outputs > 0 && ~isempty(actual_state_sequence_wf)
            state_idx_t0_wf = actual_state_sequence_wf(1);
            if state_idx_t0_wf > 0 && state_idx_t0_wf <= length(Cs) && ~isempty(Cs{state_idx_t0_wf}) && ...
               state_idx_t0_wf <= length(Ds) && ~isempty(Ds{state_idx_t0_wf})
                 C_t0 = Cs{state_idx_t0_wf}; D_t0 = Ds{state_idx_t0_wf};
                 if size(C_t0,2) == state_dim_local && size(D_t0,2) == prms.num_inputs && size(C_t0,1) == num_outputs && size(D_t0,1) == num_outputs
                     y_out(:, 1) = C_t0 * x0_ss_in + D_t0 * u_vec;
                 else
                     if ~prms.suppress_console_output; warning('generateOutputWaveforms:CDMismatchT0Dynamic', 'C/D matrix dimension mismatch for initial output with dynamic sequence.'); end
                 end
            else
                 if ~prms.suppress_console_output; warning('generateOutputWaveforms:MissingCDT0Dynamic', 'C/D model missing for initial state %d of dynamic sequence.', state_idx_t0_wf); end
            end
        end
        
        if num_points <= 1; return; end % If only t=0, we are done.

        % --- Simulate point-by-point using the dynamic sequence, repeated ---
        dynamic_seq_step_idx = 1;
        current_state_idx_wf = actual_state_sequence_wf(dynamic_seq_step_idx);
        current_interval_duration_wf = actual_interval_durations_wf(dynamic_seq_step_idx);
        time_in_current_state_wf = 0;
        time_tolerance_wf = T_period * 1e-12; % Small tolerance for time comparisons

        for k_t = 2:num_points
            time_to_propagate_this_step = t_out(k_t) - t_out(k_t-1); % dt_sim_step
            
            % Propagate through potentially multiple short dynamic intervals within one dt_sim_step
            while time_to_propagate_this_step > time_tolerance_wf
                time_left_in_current_dynamic_state = current_interval_duration_wf - time_in_current_state_wf;
                
                actual_step_duration_for_mna = 0;
                
                if time_left_in_current_dynamic_state >= time_to_propagate_this_step - time_tolerance_wf
                    % Current dt_sim_step finishes within the current dynamic interval
                    actual_step_duration_for_mna = time_to_propagate_this_step;
                    time_in_current_state_wf = time_in_current_state_wf + actual_step_duration_for_mna;
                    time_to_propagate_this_step = 0; % This dt_sim_step is done
                    % Stay in current_state_idx_wf for the next t_out point if any time remains in it
                else
                    % Current dt_sim_step is longer than what's left in the current dynamic interval
                    actual_step_duration_for_mna = time_left_in_current_dynamic_state;
                    time_in_current_state_wf = 0; % Reset for new dynamic state
                    time_to_propagate_this_step = time_to_propagate_this_step - actual_step_duration_for_mna;
                    
                    % Transition to the next dynamic state
                    dynamic_seq_step_idx = mod(dynamic_seq_step_idx, num_steps_in_dynamic_sequence) + 1;
                    current_state_idx_wf = actual_state_sequence_wf(dynamic_seq_step_idx);
                    current_interval_duration_wf = actual_interval_durations_wf(dynamic_seq_step_idx);
                end

                if actual_step_duration_for_mna > time_tolerance_wf % Only propagate if duration is significant
                    if current_state_idx_wf <= 0 || current_state_idx_wf > length(As) || isempty(As{current_state_idx_wf}) || ...
                       current_state_idx_wf > length(Bs) || isempty(Bs{current_state_idx_wf})
                         error('generateOutputWaveforms:MissingStateModelInLoop', 'State model missing for state %d at t=%.4gus.', current_state_idx_wf, t_out(k_t)*1e6);
                    end
                    A_wf = As{current_state_idx_wf}; B_wf = Bs{current_state_idx_wf};
                    
                    lambda_step_wf = calculateLambda(A_wf, actual_step_duration_for_mna, prms);
                    Phi_step_wf = eye(state_dim_local) + A_wf * lambda_step_wf;
                    if size(B_wf, 2) == prms.num_inputs
                        Gamma_step_wf = lambda_step_wf * B_wf * u_vec;
                        current_state_vec_wf = Phi_step_wf * current_state_vec_wf + Gamma_step_wf;
                    else
                        error('generateOutputWaveforms:BMismatchInLoop', 'B matrix dimension mismatch for state %d in waveform loop.', current_state_idx_wf);
                    end
                end
            end % end while time_to_propagate_this_step

            x_out(:, k_t) = current_state_vec_wf;

            if num_outputs > 0
                 if current_state_idx_wf <= 0 || current_state_idx_wf > length(Cs) || isempty(Cs{current_state_idx_wf}) || ...
                    current_state_idx_wf > length(Ds) || isempty(Ds{current_state_idx_wf})
                     if ~prms.suppress_console_output; warning('generateOutputWaveforms:MissingCDInLoop', 'C/D model missing for state %d at t=%.4gus. Output set to NaN.', current_state_idx_wf, t_out(k_t)*1e6); end
                     y_out(:, k_t) = NaN;
                 else
                     C_wf = Cs{current_state_idx_wf}; D_wf = Ds{current_state_idx_wf};
                     if size(C_wf,2) == state_dim_local && size(D_wf, 2) == prms.num_inputs && size(C_wf,1) == num_outputs && size(D_wf,1) == num_outputs
                         y_out(:, k_t) = C_wf * current_state_vec_wf + D_wf * u_vec;
                     else
                         if ~prms.suppress_console_output; warning('generateOutputWaveforms:CDMismatchInLoop', 'C/D matrix dimension mismatch for state %d. Output set to NaN.', current_state_idx_wf); end
                         y_out(:, k_t) = NaN;
                     end
                 end
            end
        end % End loop through time points k_t
        if ~prms.suppress_console_output
            fprintf('    generateOutputWaveforms: Waveform simulation finished.\n');
        end
    end % end generateOutputWaveforms



%------------------------------------------------
%------------------------------------------------
function u_vec = constructInputVector(prms)
    % Creates the input vector u based on the list of input source names
    % and values provided in config or parsed from the netlist.
    % CORRECTED: Changed prms.config_input_values to prms.input_values.

    u_var_names = prms.u_vars; % Cell array of independent source names found by parseNetlist
    num_inputs_local = length(u_var_names);
    u_vec = zeros(num_inputs_local, 1);

    % prms.input_values is the map from config (potentially empty or with overrides)
    % prms.parsed_input_values is the map of values directly parsed from the netlist
    
    config_override_values = prms.input_values; 
    netlist_parsed_values = prms.parsed_input_values;

    if num_inputs_local > 0
        if ~prms.suppress_console_output && ~isempty(keys(config_override_values))
            fprintf('    Constructing input vector u. Overrides from config.input_values will be used if present:\n');
        elseif ~prms.suppress_console_output
            fprintf('    Constructing input vector u using values parsed from netlist.\n');
        end

        for k_in = 1:num_inputs_local
            source_name = u_var_names{k_in};
            source_value = 0; % Default to 0 if not found in either map
            found_value = false;

            % Priority 1: Value from config.input_values (user override)
            if isKey(config_override_values, source_name)
                source_value = config_override_values(source_name);
                found_value = true;
                if ~prms.suppress_console_output
                    fprintf('      - Source %s: %.4g (from config.input_values override)\n', source_name, source_value);
                end
            % Priority 2: Value parsed directly from netlist
            elseif isKey(netlist_parsed_values, source_name)
                source_value = netlist_parsed_values(source_name);
                found_value = true;
                 if ~prms.suppress_console_output
                    % Only print if not overridden, to avoid duplicate messages
                    % fprintf('      - Source %s: %.4g (from netlist)\n', source_name, source_value);
                 end
            end

            if ~found_value && ~prms.suppress_console_output
                warning('InputVector:ValueNotFound', 'Value for input source "%s" not found in config.input_values or netlist parsed values. Using 0.', source_name);
            end
            u_vec(k_in) = source_value;
        end
    end
    if ~prms.suppress_console_output && num_inputs_local > 0
        fprintf('      Final u_vec = %s\n', mat2str(u_vec', 4));
    end
end % end constructInputVector


   %------------------------------------------------
  function res_struct = packageResults(prms, x0, opt_dt, dt_hist, voltage_hist, sim_time, t, x, y, valid_states_data, As, Bs, Cs, Ds, comps, sws, nodeMap_final, timing_stats_in)
        % Packages all configuration, intermediate data, and final results
        % into a single output structure for user access.
        % MODIFIED: Changed 10th parameter name to valid_states_data for clarity.

        res_struct = struct(); % Initialize the output structure

        % --- Core Results ---
        res_struct.x0_ss = x0; 
        res_struct.optimal_dt = opt_dt; 
        res_struct.dt_history = dt_hist; 
        res_struct.voltage_history = voltage_hist; 
        res_struct.simulation_time_ms = sim_time; 

        % --- Timing Results --- 
        res_struct.timing_stats = timing_stats_in; 

        % --- Waveform Data ---
        res_struct.time_vector = t; 
        res_struct.state_waveforms = x; 
        res_struct.output_waveforms = y; 
        res_struct.state_vars = prms.state_vars; 

        if iscell(prms.outputs) && isscalar(prms.outputs) && strcmp(prms.outputs{1}, 'all_states')
            res_struct.output_names = prms.state_vars; 
        else
            res_struct.output_names = prms.outputs; 
        end

        % --- Topology & Models ---
        % Line 452 (approximately, depending on exact file structure) would be here:
        res_struct.valid_states = valid_states_data; % Use the passed parameter name consistently

        if isfield(prms, 'topology_info') 
            res_struct.topology_info = prms.topology_info; 
        else
            res_struct.topology_info = struct(); 
        end
        res_struct.state_models = cell(prms.num_valid_states, 1); 
        for k_state_model = 1:prms.num_valid_states
            A_val = []; B_val = []; C_val = []; D_val = [];
            if k_state_model <= length(As) && ~isempty(As{k_state_model}); A_val = As{k_state_model}; end
            if k_state_model <= length(Bs) && ~isempty(Bs{k_state_model}); B_val = Bs{k_state_model}; end
            if k_state_model <= length(Cs) && ~isempty(Cs{k_state_model}); C_val = Cs{k_state_model}; end
            if k_state_model <= length(Ds) && ~isempty(Ds{k_state_model}); D_val = Ds{k_state_model}; end

            res_struct.state_models{k_state_model} = struct(...
                'A', A_val, ...
                'B', B_val, ...
                'C', C_val, ...
                'D', D_val);
        end

        % --- Input Configuration & Parsed Data ---
        res_struct.config = prms; 
        res_struct.components = comps; 
        res_struct.switches = sws; 
        res_struct.node_map = nodeMap_final; 

    end % end packageResults
    
   %------------------------------------------------
    function displayResults(res)
        % Displays final results summary to the console.
        % MODIFIED: Includes timing statistics if available.
        fprintf('\n========== Final Results Summary ==========\n');
        fprintf('Total Function Execution Time: %.3f ms\n', res.simulation_time_ms);

        % --- Display Timing Stats --- % <<< ADDED
        if isfield(res, 'timing_stats') && isfield(res.timing_stats, 'num_runs') && res.timing_stats.num_runs > 1
            fprintf('--- Core Calculation Timing (N=%d runs) ---\n', res.timing_stats.num_runs);
            fprintf('  Mean Time:      %.4f ms\n', res.timing_stats.mean_time_ms);
            fprintf('  Std Dev Time:   %.4f ms\n', res.timing_stats.std_dev_time_ms);
            fprintf('  Min Time:       %.4f ms\n', res.timing_stats.min_time_ms);
            fprintf('  Max Time:       %.4f ms\n', res.timing_stats.max_time_ms);
            fprintf('-------------------------------------------\n');
        end
        % --- End Timing Stats ---

        if ~isempty(res.x0_ss)
            fprintf('Steady-State Initial Conditions (x0_ss):\n');
            max_name_len = 0;
            if ~isempty(res.state_vars)
                 max_name_len = max(cellfun(@length, res.state_vars));
            end
            format_str = sprintf('  %%-%ds = %%12.6e\\n', max(max_name_len, 10)); % Dynamic padding

            for i_res = 1:length(res.x0_ss)
                 % Use a default name if state_vars is somehow shorter
                 state_name = sprintf('State %d', i_res);
                 if i_res <= length(res.state_vars) && ~isempty(res.state_vars{i_res})
                     state_name = res.state_vars{i_res};
                 end
                fprintf(format_str, state_name, res.x0_ss(i_res));
            end
        else
             fprintf('Steady-State Initial Conditions (x0_ss): N/A (No state variables)\n');
        end

        if res.config.run_optimization && ~isempty(res.optimal_dt) && ~any(isnan(res.optimal_dt))
            fprintf('Optimal Deadtimes:\n');
            dt_intervals_in_seq = 0; % Counter for deadtimes in the sequence
            seq = res.topology_info.state_sequence;
            dead_states = res.topology_info.deadtime_states;
             for k_seq_disp=1:length(seq) % Iterate through the state sequence
                 if ismember(seq(k_seq_disp), dead_states) % Check if the current state in sequence is a deadtime state
                     dt_intervals_in_seq = dt_intervals_in_seq + 1; % Increment deadtime counter
                     if dt_intervals_in_seq <= length(res.optimal_dt) % Check if we have a value for this deadtime
                         fprintf('  dt%d (State %d in seq pos %d) = %.3f ns\n', ...
                                 dt_intervals_in_seq, seq(k_seq_disp), k_seq_disp, res.optimal_dt(dt_intervals_in_seq)*1e9);
                     else
                          fprintf('  dt%d (State %d in seq pos %d) = Value Missing\n', ...
                                 dt_intervals_in_seq, seq(k_seq_disp), k_seq_disp);
                     end
                 end
             end
        elseif res.config.run_optimization
             fprintf('Optimal Deadtimes: Optimization run failed or was disabled.\n');
        else
             fprintf('Optimal Deadtimes: Optimization not run.\n');
        end
        fprintf('=========================================\n');
    end % end displayResults
    %------------------------------------------------
   function plotResults(res)
        % Plots the generated steady-state waveforms (3 periods) and
        % ZVS optimization convergence plots (if applicable).
        % MODIFIED: Adds timing performance plots if available.

        % --- Basic Checks ---
        if ~isfield(res, 'time_vector') || ~isfield(res, 'config') || isempty(res.time_vector) || (isempty(res.state_waveforms) && isempty(res.output_waveforms))
             fprintf('Skipping waveform plotting: Required data missing.\n');
             % Still check if timing plots can be made
        else
            fprintf('Plotting waveform results...\n');
            plotWaveformResults(res); % Call dedicated waveform plotting function
        end

        % --- Plot Timing Results (if available) --- % <<< ADDED
        if isfield(res, 'timing_stats') && isfield(res.timing_stats, 'run_times_ms') && ~isempty(res.timing_stats.run_times_ms)
            fprintf('Plotting timing results...\n');
            plotTimingResults(res.timing_stats); % Call dedicated timing plotting function
        else
             fprintf('Skipping timing plots: No timing data available.\n');
        end

    end % end plotResults
    % --- Nested Function for Waveform Plots ---
    % --- Nested Function for Waveform Plots ---
    function plotWaveformResults(res)
        % Plots the generated steady-state waveforms (3 periods) and
        % ZVS optimization convergence plots (if applicable).
        % MODIFIED: Uses generateSequenceFromTimings for accurate switch state plotting.
        % MODIFIED: ZVS convergence plots now reflect per-switch deadtimes/targets.

        % --- Basic Checks ---
        if ~isfield(res, 'time_vector') || ~isfield(res, 'config') || isempty(res.time_vector) 
            fprintf('Skipping waveform plotting: Required data (time_vector/config) missing.\n');
            return;
        end
        if isempty(res.state_waveforms) && isempty(res.output_waveforms) && res.config.num_switches == 0
            fprintf('Skipping waveform plotting: No state/output waveforms and no switches.\n');
            return;
        end
        
        fprintf('Plotting waveform results...\n');

        % --- Determine if Convergence Plots are Needed ---
        plot_convergence = false;
        num_switches_for_convergence_plot = 0;

        if isfield(res.config, 'run_optimization') && res.config.run_optimization && ...
           isfield(res, 'dt_history') && ~isempty(res.dt_history) && size(res.dt_history, 1) > 0 && ...
           isfield(res, 'voltage_history') && ~isempty(res.voltage_history) && size(res.voltage_history, 1) > 0
            
            if size(res.dt_history, 2) == res.config.num_switches && ...
               size(res.voltage_history, 2) == res.config.num_switches && ...
               res.config.num_switches > 0
                plot_convergence = true;
                num_switches_for_convergence_plot = res.config.num_switches;
                fprintf('  (Including ZVS optimization convergence plots for %d switches)\n', num_switches_for_convergence_plot);
            else
                warning('PlotWaveformResults:ConvergenceHistoryMismatch', ...
                        'Size of dt_history/voltage_history columns (%d/%d) does not match number of switches (%d). Skipping ZVS convergence plots.', ...
                        size(res.dt_history, 2), size(res.voltage_history, 2), res.config.num_switches);
                plot_convergence = false;
            end
        end

        % --- Setup Figure and Subplots ---
        num_state_plots = 0;
        if isfield(res, 'state_waveforms') && ~isempty(res.state_waveforms)
            num_state_plots = size(res.state_waveforms, 1);
        end
        
        num_output_plots = 0;
        if isfield(res, 'output_waveforms') && ~isempty(res.output_waveforms)
            num_output_plots = size(res.output_waveforms, 1);
        end

        num_plots_wave = 1; % Start with 1 for switch states (even if no switches, we'll show a placeholder)

        is_all_states_only = false;
        if iscell(res.config.outputs) && isscalar(res.config.outputs) && strcmp(res.config.outputs{1}, 'all_states')
            is_all_states_only = true;
            if num_state_plots > 0
                num_plots_wave = 1 + num_state_plots; 
            end
        elseif num_output_plots > 0
             num_plots_wave = 1 + num_output_plots;
        end
        
        num_plots_total = num_plots_wave;
        if plot_convergence
            num_plots_total = num_plots_total + 2; % Add 2 plots for convergence
        end
        
        if num_plots_total <= 1 && res.config.num_switches > 0 % If only switch plot is planned
             num_plots_total = 1;
        elseif num_plots_total == 0 % No plots at all
            fprintf('  No data to plot for waveforms or convergence.\n');
            return;
        end


        fig_handle_wave = figure;
        set(fig_handle_wave, 'Name', 'AnalyzeSwitchingConverter Waveform Results', 'NumberTitle', 'off');
        ax_wave = gobjects(max(1,num_plots_total), 1); % Ensure at least 1 axis

        % --- Plotting Parameters ---
        t_us = res.time_vector * 1e6;
        T_us = res.config.T_period * 1e6;
        num_sw = res.config.num_switches;
        plot_idx = 1;

        num_colors_needed = max([num_sw, num_state_plots, num_output_plots, ifelse(plot_convergence, num_switches_for_convergence_plot, 0)]);
        if num_colors_needed == 0; num_colors_needed = 1; end
        plot_colors = lines(max(1,num_colors_needed)); % Ensure lines() gets at least 1

        try 
            % --- Plot 1: Switching States ---
            ax_wave(plot_idx) = subplot(num_plots_total, 1, plot_idx);
            hold(ax_wave(plot_idx), 'on'); grid(ax_wave(plot_idx), 'on'); box(ax_wave(plot_idx), 'on');
            title(ax_wave(plot_idx), 'Switching States (Event-Driven Reconstruction)');
            ylabel(ax_wave(plot_idx), 'State');

            if num_sw > 0 && isfield(res, 'valid_states') && ~isempty(res.valid_states)
                % Reconstruct sequence and durations using final deadtimes
                prms_plot = res.config; % Start with the full config
                prms_plot.deadtimes_applied = res.optimal_dt; % Use final per-switch deadtimes
                if isempty(prms_plot.deadtimes_applied) || any(isnan(prms_plot.deadtimes_applied))
                    % Fallback if optimal_dt is not good
                    prms_plot.deadtimes_applied = res.config.initial_deadtimes; 
                    warning('PlotWaveformResults:UsingInitialDTForPlot', 'Optimal deadtimes invalid for plotting, using initial_deadtimes for switch state reconstruction.');
                end
                if length(prms_plot.deadtimes_applied) ~= num_sw
                     warning('PlotWaveformResults:PlotDTLengthMismatch', 'Length of deadtimes for plotting (%d) does not match num_switches (%d). Plot may be inaccurate.', length(prms_plot.deadtimes_applied), num_sw);
                     % Attempt to fix if it's a scalar and num_sw > 0
                     if isscalar(prms_plot.deadtimes_applied) && num_sw > 0
                         prms_plot.deadtimes_applied = repmat(prms_plot.deadtimes_applied, 1, num_sw);
                     else
                         prms_plot.deadtimes_applied = zeros(1, num_sw); % Fallback to zero
                     end
                end

                prms_plot.suppress_console_output = true; % Suppress messages from generateSequenceFromTimings during plotting

                [plot_state_sequence, plot_interval_durations] = generateSequenceFromTimings(prms_plot);
                
                valid_states_matrix = res.valid_states;
                switch_waveforms = zeros(num_sw, length(t_us));
                time_tolerance_wf_plot = T_us * 1e-9; 

                if ~isempty(plot_state_sequence)
                    current_dynamic_step_plot = 1;
                    current_dynamic_state_idx_plot = plot_state_sequence(current_dynamic_step_plot);
                    current_dynamic_duration_us_plot = plot_interval_durations(current_dynamic_step_plot) * 1e6;
                    time_in_current_dynamic_state_us_plot = 0;

                    for k_t_plt = 1:length(t_us)
                        t_point_us = t_us(k_t_plt);
                        
                        % Determine which dynamic interval this t_point falls into
                        % This logic assumes t_us is monotonically increasing and covers multiple periods
                        time_within_period_us = mod(t_point_us, T_us);
                        % Handle endpoint case where mod(T_us, T_us) is 0 but should be T_us
                        if abs(time_within_period_us) < time_tolerance_wf_plot && t_point_us > time_tolerance_wf_plot 
                            if abs(t_point_us - round(t_point_us/T_us)*T_us) < time_tolerance_wf_plot
                                time_within_period_us = T_us;
                            end
                        end

                        cumulative_duration_us = 0;
                        active_plot_state_for_t = 0;
                        for k_seq_plot = 1:length(plot_state_sequence)
                            duration_this_interval_us = plot_interval_durations(k_seq_plot) * 1e6;
                            if time_within_period_us < cumulative_duration_us + duration_this_interval_us + time_tolerance_wf_plot
                                active_plot_state_for_t = plot_state_sequence(k_seq_plot);
                                break;
                            end
                            cumulative_duration_us = cumulative_duration_us + duration_this_interval_us;
                        end
                        if active_plot_state_for_t == 0 % Should fall into the last state if at T_us
                             active_plot_state_for_t = plot_state_sequence(end);
                        end
                        
                        if active_plot_state_for_t > 0 && active_plot_state_for_t <= size(valid_states_matrix, 1)
                            switch_waveforms(:, k_t_plt) = valid_states_matrix(active_plot_state_for_t, :)';
                        else
                            switch_waveforms(:, k_t_plt) = zeros(num_sw, 1); 
                        end
                    end

                    stairs_handles = gobjects(1, num_sw); legend_entries = cell(1, num_sw);
                    for sw_plt = 1:num_sw
                        stairs_handles(sw_plt) = stairs(ax_wave(plot_idx), t_us, switch_waveforms(sw_plt, :) + (sw_plt-1)*1.1, 'Color', plot_colors(sw_plt,:), 'LineWidth', 1.5);
                        legend_entries{sw_plt} = res.config.switches_ordered_list(sw_plt).name;
                    end
                    ylim(ax_wave(plot_idx), [-0.1, max(1, num_sw*1.1 + 0.1)]);
                    if num_sw > 0
                        yticks(ax_wave(plot_idx), 0.5 + (0:(num_sw-1))*1.1);
                        yticklabels(ax_wave(plot_idx), {res.config.switches_ordered_list.name});
                        legend(stairs_handles, legend_entries, 'Location', 'bestoutside');
                    end
                else
                    title(ax_wave(plot_idx), 'Switching States (Failed to reconstruct sequence for plot)');
                end
            else
                title(ax_wave(plot_idx), 'Switching States (No Switches in Netlist)');
                text(ax_wave(plot_idx), mean(get(ax_wave(plot_idx),'XLim')), mean(get(ax_wave(plot_idx),'YLim')), 'No Switches in Circuit', 'HorizontalAlignment', 'center');
            end
            plot_idx = plot_idx + 1;

            % --- Plot State Variables or Outputs ---
            if is_all_states_only && num_state_plots > 0
                for i_state_plot = 1:num_state_plots
                    ax_wave(plot_idx) = subplot(num_plots_total, 1, plot_idx);
                    state_name = res.state_vars{i_state_plot};
                    plot_color = plot_colors(mod(i_state_plot-1, size(plot_colors,1))+1,:);
                    plot(ax_wave(plot_idx), t_us, res.state_waveforms(i_state_plot, :), 'Color', plot_color, 'LineWidth', 1.5);
                    grid(ax_wave(plot_idx), 'on'); box(ax_wave(plot_idx), 'on');
                    ylabel(ax_wave(plot_idx), strrep(state_name, '_', '\_')); 
                    plot_idx = plot_idx + 1;
                end
            elseif num_output_plots > 0 && ~is_all_states_only
                for i_out_plot = 1:num_output_plots
                     ax_wave(plot_idx) = subplot(num_plots_total, 1, plot_idx);
                     out_name = res.output_names{i_out_plot};
                     plot_color = plot_colors(mod(i_out_plot-1, size(plot_colors,1))+1,:); 
                     plot(ax_wave(plot_idx), t_us, res.output_waveforms(i_out_plot, :), 'Color', plot_color, 'LineWidth', 1.5);
                     grid(ax_wave(plot_idx), 'on'); box(ax_wave(plot_idx), 'on');
                     ylabel(ax_wave(plot_idx), strrep(out_name, '_', '\_')); 
                     plot_idx = plot_idx + 1;
                end
            end

            % --- Plot Convergence Plots (Conditional) ---
            if plot_convergence && num_switches_for_convergence_plot > 0
                num_iters = size(res.dt_history, 1);
                iters = 1:num_iters;

                if num_iters > 0
                    % Deadtime Convergence Plot
                    ax_wave(plot_idx) = subplot(num_plots_total, 1, plot_idx);
                    hold(ax_wave(plot_idx), 'on'); grid(ax_wave(plot_idx), 'on'); box(ax_wave(plot_idx), 'on');
                    title(ax_wave(plot_idx), 'ZVS Deadtime Optimization Convergence (Per Switch)');
                    ylabel(ax_wave(plot_idx), 'Deadtime (ns)');
                    dt_handles = gobjects(1, num_switches_for_convergence_plot); 
                    dt_legends = cell(1, num_switches_for_convergence_plot);
                    for i_sw_conv = 1:num_switches_for_convergence_plot
                        color_idx = mod(i_sw_conv-1, size(plot_colors,1)) + 1;
                        dt_handles(i_sw_conv) = plot(ax_wave(plot_idx), iters, res.dt_history(:, i_sw_conv) * 1e9, '.-', 'MarkerSize', 10, 'Color', plot_colors(color_idx,:));
                        dt_legends{i_sw_conv} = sprintf('dt for %s', res.config.switches_ordered_list(i_sw_conv).name);
                    end
                    if num_iters == 1; xticks(ax_wave(plot_idx), 1); end
                    legend(dt_handles, dt_legends, 'Location', 'best');
                    plot_idx = plot_idx + 1;

                    % Target Voltage Convergence Plot
                    ax_wave(plot_idx) = subplot(num_plots_total, 1, plot_idx);
                    hold(ax_wave(plot_idx), 'on'); grid(ax_wave(plot_idx), 'on'); box(ax_wave(plot_idx), 'on');
                    title(ax_wave(plot_idx), 'ZVS Target Voltage Convergence (Per Targeted Switch)');
                    ylabel(ax_wave(plot_idx), 'Target Voltage (V)');
                    xlabel(ax_wave(plot_idx), 'Iteration'); 
                    vc_handles = gobjects(1, num_switches_for_convergence_plot); 
                    vc_legends = cell(1, num_switches_for_convergence_plot);
                    
                    actual_targets_plotted = 0;
                    for i_sw_conv = 1:num_switches_for_convergence_plot
                        % Check if this switch was actually targeted and has a valid cap index
                        is_targeted = false;
                        cap_var_name_for_legend = 'N/A';
                        if isfield(res.config,'is_switch_targeted_for_zvs') && length(res.config.is_switch_targeted_for_zvs) >= i_sw_conv && res.config.is_switch_targeted_for_zvs(i_sw_conv)
                           is_targeted = true;
                           if isfield(res.config,'zvs_cap_state_var_indices') && length(res.config.zvs_cap_state_var_indices) >= i_sw_conv
                               cap_idx_for_legend = res.config.zvs_cap_state_var_indices(i_sw_conv);
                               if ~isnan(cap_idx_for_legend) && cap_idx_for_legend > 0 && cap_idx_for_legend <= length(res.state_vars)
                                   cap_var_name_for_legend = res.state_vars{cap_idx_for_legend};
                               end
                           end
                        end
                        
                        if ~all(isnan(res.voltage_history(:, i_sw_conv))) % Only plot if there's data
                            actual_targets_plotted = actual_targets_plotted + 1;
                            color_idx = mod(i_sw_conv-1, size(plot_colors,1)) + 1;
                            vc_handles(actual_targets_plotted) = plot(ax_wave(plot_idx), iters, res.voltage_history(:, i_sw_conv), '.-', 'MarkerSize', 10, 'Color', plot_colors(color_idx,:));
                            vc_legends{actual_targets_plotted} = sprintf('Vc for %s (%s)', res.config.switches_ordered_list(i_sw_conv).name, strrep(cap_var_name_for_legend,'_','\_'));
                        end
                    end
                    yline(ax_wave(plot_idx), res.config.voltage_tolerance, '--r', 'HandleVisibility', 'off');
                    yline(ax_wave(plot_idx), -res.config.voltage_tolerance, '--r', 'HandleVisibility', 'off');
                    if num_iters == 1; xticks(ax_wave(plot_idx), 1); end
                    if actual_targets_plotted > 0
                        legend(vc_handles(1:actual_targets_plotted), vc_legends(1:actual_targets_plotted), 'Location', 'best');
                    end
                    plot_idx = plot_idx + 1;
                end
            end 

            % --- Final Figure Adjustments ---
            for i_ax = 1:(plot_idx - 2) % All subplots except the last one plotted
                if isgraphics(ax_wave(i_ax)); set(ax_wave(i_ax), 'XTickLabel', []); end
            end
            if plot_idx > 1 && isgraphics(ax_wave(plot_idx-1))
                 if plot_convergence && (plot_idx-1 == num_plots_wave + 2) 
                     % Last plot was voltage convergence, xlabel already set
                 else 
                     xlabel(ax_wave(plot_idx-1), 'Time (\mus)');
                 end
            end

            link_indices = 1:max(1, num_plots_wave); % Link only waveform plots
            valid_axes_to_link = ax_wave(link_indices(isgraphics(ax_wave(link_indices))));
            if ~isempty(valid_axes_to_link)
                linkaxes(valid_axes_to_link, 'x');
                if ~isempty(t_us); xlim(valid_axes_to_link(1), [0 t_us(end)]); end
            end

            sgtitle(fig_handle_wave, sprintf('Steady-State Analysis Waveforms (fs=%.1fkHz)', res.config.fs/1e3), 'FontWeight', 'bold');

        catch ME_plot
             warning('AnalyzeSwitchingConverter:PlottingError', 'Error during waveform plotting: %s\nIn %s (line %d)', ME_plot.message, ME_plot.stack(1).name, ME_plot.stack(1).line);
             if exist('fig_handle_wave', 'var') && ishandle(fig_handle_wave); close(fig_handle_wave); end
        end
    end % end plotWaveformResults

    % --- NEW Nested Function for Timing Plots ---
    function plotTimingResults(timing_stats)
        % Creates a new figure and plots the timing results:
        % 1. Simulation time vs. run number
        % 2. Histogram of simulation times

        if isempty(timing_stats) || ~isfield(timing_stats, 'run_times_ms') || isempty(timing_stats.run_times_ms) || timing_stats.num_runs <= 1
            fprintf('  Insufficient data for timing plots.\n');
            return;
        end

        run_times = timing_stats.run_times_ms;
        num_runs = timing_stats.num_runs;
        mean_time = timing_stats.mean_time_ms;
        run_indices = 1:num_runs;

        fig_handle_time = figure; % Create a new figure for timing
        set(fig_handle_time, 'Name', 'AnalyzeSwitchingConverter Timing Results', 'NumberTitle', 'off');

        % Subplot 1: Time vs Run Number
        ax1 = subplot(2, 1, 1);
        plot(ax1, run_indices, run_times, 'b.-', 'MarkerSize', 8);
        hold(ax1, 'on');
        plot(ax1, run_indices, repmat(mean_time, 1, num_runs), 'r--', 'LineWidth', 1.5);
        hold(ax1, 'off');
        grid(ax1, 'on');
        xlabel(ax1, 'Run number');
        ylabel(ax1, 'Simulation time (ms)');
        title(ax1, sprintf('ESVA Simulation Times (%d runs)', num_runs));
        legend(ax1, 'Individual runs', sprintf('Average (%.2f ms)', mean_time), 'Location', 'best');
        if num_runs > 1
            xlim(ax1, [0 num_runs + 1]); % Adjust xlim for better visibility
        end
        if ~isempty(run_times) % Prevent error if run_times is empty
           ylim(ax1, [0, max(run_times)*1.1 + eps]); % Adjust ylim slightly above max
        end


        % Subplot 2: Histogram
        ax2 = subplot(2, 1, 2);
        histogram(ax2, run_times);
        hold(ax2, 'on');
        % Add vertical line for mean
        plot(ax2, [mean_time mean_time], get(ax2,'YLim'), 'r--', 'LineWidth', 1.5, 'HandleVisibility','off');
        hold(ax2, 'off');
        grid(ax2, 'on');
        xlabel(ax2, 'Simulation time (ms)');
        ylabel(ax2, 'Frequency');
        title(ax2, 'Distribution of Simulation Times');

    end % end plotTimingResults
    % --- End New Plotting Functions ---

    %------------------------------------------------
    function result = ifelse(condition, true_val, false_val)
        % Helper function for conditional assignment (like ternary operator)
        if condition
            result = true_val;
        else
            result = false_val;
        end
    end
    %------------------------------------------------

    %% --- Event-Driven Sequence Generation Method ---
       function [actual_state_sequence, actual_interval_durations] = generateSequenceFromTimings(prms)
        % generateSequenceFromTimings: Generates the state sequence and durations 
        % based on an event-driven timing method.
        % This method processes ideal per-switch ON/OFF timings and applies
        % per-switch deadtimes to determine the actual sequence of operations
        % and their precise durations over one switching period.
        % This function is designed to be N-switch adaptable.
        %
        % Inputs:
        %   prms - Struct containing all necessary parameters:
        %     .ts_switches         - Cell array: {i} = vector of ideal ON/OFF instants for switch i.
        %                            Assumed format: [on1, off1, on2, off2, ...]. Length must match num_switches.
        %     .deadtimes_applied   - Vector: deadtimes_applied(i) = deadtime for switch i's ON transitions. Length must match num_switches.
        %     .num_switches        - Scalar: Total number of switches (N).
        %     .T_period            - Scalar: Switching period (1/fs).
        %     .valid_states        - Matrix: Rows are valid switch state vectors.
        %                            Columns correspond to switches in prms.switches_ordered_list order. Number of columns must match num_switches.
        %     .switches_ordered_list - Struct array: Info about switches, including 'name'. Order defines switch correspondence. Length must match num_switches.
        %     .time_tolerance_uniquetol - Scalar: Absolute tolerance for uniquetol().
        %     .min_interval_duration    - Scalar: Minimum duration for an interval to be considered significant.
        %     .suppress_console_output - Logical: Flag to suppress fprintf messages.
        %
        % Outputs:
        %   actual_state_sequence     - Row vector: Sequence of state indices (from prms.valid_states).
        %   actual_interval_durations - Row vector: Durations of each state in the sequence.

        if ~prms.suppress_console_output
            fprintf('    Generating sequence and durations using event-driven timing method...\n');
        end
        
        num_total_switches = prms.num_switches; % N switches
        T_period_in = prms.T_period;
        ts_switches_ideal_in = prms.ts_switches;
        deadtimes_applied_in = prms.deadtimes_applied;
        valid_states_matrix_in = prms.valid_states;
        min_interval_tol_in = prms.min_interval_duration;
        time_tol_uniquetol_abs = prms.time_tolerance_uniquetol;

        actual_state_sequence = [];
        actual_interval_durations = [];
        
        % Handle case with no switches
        if num_total_switches == 0
            if isempty(valid_states_matrix_in) || (size(valid_states_matrix_in,1)==1 && all(valid_states_matrix_in == 0))
                 actual_state_sequence = 1; 
            else
                 all_off_row = find(all(valid_states_matrix_in == 0, 2), 1);
                 if ~isempty(all_off_row)
                     actual_state_sequence = all_off_row;
                 else
                     actual_state_sequence = 1; 
                 end
            end
            actual_interval_durations = T_period_in;
            if ~prms.suppress_console_output
                fprintf('    No switches. Sequence is static state %d for full period.\n', actual_state_sequence);
            end
            return;
        end

        % --- Input Validation (already performed more extensively in main function body) ---
        % Basic checks here are for robustness if called directly with a malformed prms.
        if ~iscell(ts_switches_ideal_in) || length(ts_switches_ideal_in) ~= num_total_switches
            error('generateSequenceFromTimings:InputMismatch_ts_ideal', 'prms.ts_switches must be a cell array of length %d (number of switches).', num_total_switches);
        end
        if ~isnumeric(deadtimes_applied_in) || ~isvector(deadtimes_applied_in) || length(deadtimes_applied_in) ~= num_total_switches
            error('generateSequenceFromTimings:InputMismatch_dt_applied', 'prms.deadtimes_applied must be a numeric vector of length %d (number of switches).', num_total_switches);
        end
        if any(deadtimes_applied_in < 0)
            error('generateSequenceFromTimings:NegativeDeadtime', 'Applied deadtimes must be non-negative.');
        end
        if T_period_in <= 0
            error('generateSequenceFromTimings:InvalidPeriod', 'Switching period T_period must be positive.');
        end
        if isempty(valid_states_matrix_in) && num_total_switches > 0 % Allow empty if no switches, handled above
            error('generateSequenceFromTimings:MissingValidStates', 'prms.valid_states matrix is empty but switches exist.');
        end
        if ~isempty(valid_states_matrix_in) && size(valid_states_matrix_in, 2) ~= num_total_switches
            error('generateSequenceFromTimings:ValidStatesColMismatch', 'Number of columns in prms.valid_states (%d) must match number of switches (%d).', size(valid_states_matrix_in, 2), num_total_switches);
        end
        if length(prms.switches_ordered_list) ~= num_total_switches
             error('generateSequenceFromTimings:SwitchesOrderedListMismatch', 'Length of prms.switches_ordered_list (%d) must match num_switches (%d).', length(prms.switches_ordered_list), num_total_switches);
        end


        % --- a. Adjust Ideal Timings for Deadtime ---
        t_prime_switches_adjusted = cell(1, num_total_switches);
        all_event_times_list = cell(1, num_total_switches); 

        for i_sw = 1:num_total_switches % Loop for N switches
            tsi = ts_switches_ideal_in{i_sw}; 
            dti = deadtimes_applied_in(i_sw); 
            t_prime_si_current = tsi;       

            if isempty(tsi) 
                t_prime_switches_adjusted{i_sw} = [];
                all_event_times_list{i_sw} = []; 
                continue; 
            end
            
            if ~issorted(tsi) || any(tsi < -time_tol_uniquetol_abs) || any(tsi > T_period_in + time_tol_uniquetol_abs)
                if ~prms.suppress_console_output
                    warning('generateSequenceFromTimings:UnsortedOrOutOfBoundIdealTimes', ...
                            'Ideal times for switch %d ("%s") not sorted or out of [0, T_period]. Clamping/Sorting. Times: [%s]', ...
                            i_sw, prms.switches_ordered_list(i_sw).name, num2str(tsi));
                end
                tsi = sort(max(0, min(tsi, T_period_in))); 
                t_prime_si_current = tsi; 
            end
            
            for j_edge = 1:length(tsi)
                if mod(j_edge, 2) == 1 
                    t_prime_si_current(j_edge) = tsi(j_edge) + dti;
                end
            end
            t_prime_switches_adjusted{i_sw} = t_prime_si_current;
            all_event_times_list{i_sw} = t_prime_si_current; 
        end
        
        all_adjusted_event_times_flat = horzcat(all_event_times_list{:}); 

        % --- b. Merge and Sort All Event Times into a Master Timeline ---
        master_event_timeline = [0, all_adjusted_event_times_flat, T_period_in];
        
        master_event_timeline(master_event_timeline < 0 & abs(master_event_timeline) > time_tol_uniquetol_abs) = 0; 
        master_event_timeline(master_event_timeline < 0) = 0; 
        master_event_timeline(master_event_timeline > T_period_in & (master_event_timeline - T_period_in) > time_tol_uniquetol_abs) = T_period_in; 
        master_event_timeline(master_event_timeline > T_period_in) = T_period_in; 

        ts_merged_unique = uniquetol(sort(master_event_timeline), time_tol_uniquetol_abs);

        if isempty(ts_merged_unique) || abs(ts_merged_unique(1) - 0) > time_tol_uniquetol_abs / 2
            ts_merged_unique = [0; ts_merged_unique(:)]; 
        else
            ts_merged_unique(1) = 0; 
        end
        if abs(ts_merged_unique(end) - T_period_in) > time_tol_uniquetol_abs / 2
            if ts_merged_unique(end) < T_period_in - time_tol_uniquetol_abs / 2 
                ts_merged_unique = [ts_merged_unique(:); T_period_in]; 
            else 
                 ts_merged_unique(end) = T_period_in; 
            end
        else
            ts_merged_unique(end) = T_period_in; 
        end
        ts_merged_unique = uniquetol(sort(ts_merged_unique), time_tol_uniquetol_abs);
        ts_merged_unique = ts_merged_unique(ts_merged_unique >= 0 & ts_merged_unique <= T_period_in + time_tol_uniquetol_abs); 
        ts_merged_unique(ts_merged_unique > T_period_in) = T_period_in; 
        ts_merged_unique = uniquetol(ts_merged_unique, time_tol_uniquetol_abs); 


        % --- c. Determine Switch States for Each Interval in the Master Timeline ---
        actual_state_sequence_temp = [];
        actual_interval_durations_temp = [];

        if length(ts_merged_unique) < 2 
             if ~prms.suppress_console_output
                fprintf('    generateSequenceFromTimings: Merged timeline has < 2 points. Evaluating state at t=T/2 for a single interval.\n');
             end
            t_eval_single = T_period_in / 2; 
            current_switch_vector_single = zeros(1, num_total_switches); % Sized for N switches
            for i_sw = 1:num_total_switches % Loop for N switches
                t_prime_si_adj = t_prime_switches_adjusted{i_sw}; 
                current_switch_vector_single(i_sw) = 0; 
                if ~isempty(t_prime_si_adj)
                    for j_edge = 1:2:length(t_prime_si_adj) 
                        on_time_adj = t_prime_si_adj(j_edge);
                        off_time_adj = T_period_in + time_tol_uniquetol_abs; 
                        if j_edge + 1 <= length(t_prime_si_adj)
                            off_time_adj = t_prime_si_adj(j_edge+1);
                        end
                        if t_eval_single >= on_time_adj && t_eval_single < off_time_adj
                            current_switch_vector_single(i_sw) = 1; 
                            break; 
                        end
                    end
                end
            end
            [is_match_single, state_idx_single] = ismember(current_switch_vector_single, valid_states_matrix_in, 'rows');
            if is_match_single
                actual_state_sequence_temp = state_idx_single;
                actual_interval_durations_temp = T_period_in;
            else
                error('generateSequenceFromTimings:InvalidStateSingleInt', 'Generated switch vector [%s] for single interval not found in valid_states_matrix_in.', num2str(current_switch_vector_single));
            end
        else 
            for k_interval = 1:(length(ts_merged_unique) - 1)
                t_start = ts_merged_unique(k_interval);
                t_end = ts_merged_unique(k_interval + 1);
                duration_k = t_end - t_start;

                if duration_k < min_interval_tol_in && duration_k > eps(T_period_in) 
                    if ~prms.suppress_console_output
                        fprintf('    Skipping interval (%.6es to %.6es, duration: %.3es) as it''s < min_interval_duration (%.3es).\n', t_start, t_end, duration_k, min_interval_tol_in);
                    end
                    continue; 
                elseif duration_k <= eps(T_period_in) 
                    continue;
                end

                t_eval = (t_start + t_end) / 2; 
                current_switch_vector = zeros(1, num_total_switches); % Sized for N switches

                for i_sw = 1:num_total_switches % Loop for N switches
                    t_prime_si_adj = t_prime_switches_adjusted{i_sw}; 
                    current_switch_vector(i_sw) = 0; 
                    if ~isempty(t_prime_si_adj)
                        for j_edge = 1:2:length(t_prime_si_adj) 
                            on_time_adj = t_prime_si_adj(j_edge);
                            off_time_adj = T_period_in + time_tol_uniquetol_abs; 
                            if j_edge + 1 <= length(t_prime_si_adj)
                                off_time_adj = t_prime_si_adj(j_edge+1);
                            end
                            if t_eval >= on_time_adj && t_eval < off_time_adj
                                current_switch_vector(i_sw) = 1; 
                                break; 
                            end
                        end
                    end
                end
                
                [is_match, state_idx_k] = ismember(current_switch_vector, valid_states_matrix_in, 'rows');
                if ~is_match
                    error('generateSequenceFromTimings:InvalidStateGenerated', ...
                          'Generated switch vector [%s] at t_eval=%.6es (interval [%.6es, %.6es]) does not map to any row in valid_states_matrix_in. Check ideal timings, deadtimes, and valid_states definitions.', ...
                          num2str(current_switch_vector), t_eval, t_start, t_end);
                end
                
                actual_state_sequence_temp = [actual_state_sequence_temp, state_idx_k];
                actual_interval_durations_temp = [actual_interval_durations_temp, duration_k];
            end
        end
        
        % --- d. Consolidate Consecutive Identical States in the Sequence ---
        if length(actual_state_sequence_temp) > 1
            actual_state_sequence = actual_state_sequence_temp(1);
            actual_interval_durations = actual_interval_durations_temp(1);
            for i_seq = 2:length(actual_state_sequence_temp)
                if actual_state_sequence_temp(i_seq) == actual_state_sequence(end)
                    actual_interval_durations(end) = actual_interval_durations(end) + actual_interval_durations_temp(i_seq);
                else
                    actual_state_sequence = [actual_state_sequence, actual_state_sequence_temp(i_seq)];
                    actual_interval_durations = [actual_interval_durations, actual_interval_durations_temp(i_seq)];
                end
            end
        else 
            actual_state_sequence = actual_state_sequence_temp;
            actual_interval_durations = actual_interval_durations_temp;
        end

        % --- Final Sanity Checks and Normalization ---
        if isempty(actual_state_sequence) && T_period_in > 0
            if ~prms.suppress_console_output
                warning('generateSequenceFromTimings:EmptySequencePostConsolidation', 'Sequence is empty after consolidation. This might happen if all generated intervals are shorter than min_interval_duration or deadtimes are too large. Defaulting to an all-off state for the full period.');
            end
            all_off_vector = zeros(1, num_total_switches); % Sized for N switches
            [is_member_off, state_idx_off] = ismember(all_off_vector, valid_states_matrix_in, 'rows');
            if is_member_off
                actual_state_sequence = state_idx_off;
            else
                actual_state_sequence = 1; 
                if ~prms.suppress_console_output && num_total_switches > 0
                     warning('generateSequenceFromTimings:NoAllOffStateFoundForFallback', 'All-off state not found in valid_states for fallback. Using state 1.');
                end
            end
            actual_interval_durations = T_period_in;
        end

        current_total_duration = sum(actual_interval_durations);
        if abs(current_total_duration - T_period_in) > time_tol_uniquetol_abs && current_total_duration > time_tol_uniquetol_abs 
            if ~prms.suppress_console_output
                fprintf('    Normalizing total duration of generated sequence from %.6es to %.6es.\n', current_total_duration, T_period_in);
            end
            actual_interval_durations = actual_interval_durations * (T_period_in / current_total_duration);
        end
        actual_interval_durations(actual_interval_durations < 0) = 0;

        if ~prms.suppress_console_output
            fprintf('    Event-driven method generated sequence: [%s]\n', num2str(actual_state_sequence));
            fprintf('    Corresponding durations (us): [%s]\n', sprintf('%.4f ', actual_interval_durations*1e6));
            fprintf('    Sum of generated durations: %.6f us (Target: %.6f us)\n', sum(actual_interval_durations)*1e6, T_period_in*1e6);
        end

    end % end generateSequenceFromTimings



end % end AnalyzeSwitchingConverter (main function)
