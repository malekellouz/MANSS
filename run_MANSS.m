% =========================================================================
% Runner Script for AnalyzeSwitchingConverter (Cuk Converter Version)
% =========================================================================
% This script is adapted from the buck converter runner to analyze the
% synchronous Cuk converter.
%
% Workflow:
% 1. Save the Cuk Converter netlist text into a file (e.g., 'Cuk_netlist.txt').
% 2. Run Section 1 to parse the netlist and identify the switches (M1, M2).
% 3. Run Section 2 to configure the specific timings for the Cuk converter.
% 4. Run Section 3 to execute the full steady-state analysis.
% =========================================================================
clear; clc; close all;
fprintf('--- Starting Analysis Runner Script (Cuk Converter Version) ---\n');

%% =========================================================================
% SECTION 1: Define Common Parameters, Parse Netlist & Identify Switches
% =========================================================================
fprintf('\n--- Running Section 1: Define Common Params, Parse Netlist & Identify Switches ---\n');

% --- Define Netlist File Path ---
% Point to the netlist file.
netlistFilePath = 'cuk_netlist.txt';
fprintf('Using netlist file: %s\n', netlistFilePath);
if ~isfile(netlistFilePath); error('RunnerScript:NetlistNotFound', 'Netlist file not found: %s. Please save the netlist content.', netlistFilePath); end
% --- Define Essential Physical Parameters ---
common_params = struct();
common_params.Ron = 0.1;  % Default ON resistance for switches
common_params.Roff = 1e6; % Default OFF resistance for switches
common_params.fs = 25e3; % Switching frequency (Hz)
common_params.T_period = 1/common_params.fs; % Period (s)

% --- Define Input Source Values (Optional Overrides) ---
% If specific values for V or I sources in the netlist are desired, define them here.
% Otherwise, values from the netlist will be used.
% Format: common_params.input_values('SourceName') = Value;
common_params.input_values = containers.Map('KeyType','char', 'ValueType','double');
% Example for buck_netlist_2.txt:
% common_params.input_values('V1') = 100.0; % Overrides V1 value if it's different in netlist

% --- Define Parasitic Modeling Option ---
common_params.include_parasitics = false; % Set to true to include Rs and Cs for switches

% --- Define Individual Switch Parasitics (Required if include_parasitics=true) ---
% USER ACTION REQUIRED (potentially after first run of Section 1):
% After running Section 1 once, 'switches_identified_ordered' will show the
% names of switches from YOUR netlist (e.g., M1, M2, S_A, etc.).
% You MUST then define parasitics for EACH of these identified switches here.
% The field names (e.g., .M1, .S_A) MUST EXACTLY MATCH the switch names.
common_params.switch_parasitics = struct();
if common_params.include_parasitics
    fprintf('INFO: Parasitics are included. Ensure common_params.switch_parasitics is correctly populated for switches identified from YOUR netlist.\n');
    % Example structure (REPLACE WITH YOUR SWITCH NAMES AND VALUES):
    % if isfield(switches_identified_ordered(1), 'name') % Check if first switch exists
    %    switch_name_1 = switches_identified_ordered(1).name; % e.g., 'M1'
    %    common_params.switch_parasitics.(switch_name_1) = struct('Cs', 2e-9, 'Rs', 0.2);
    % end
    % if length(switches_identified_ordered) >= 2 && isfield(switches_identified_ordered(2), 'name')
    %    switch_name_2 = switches_identified_ordered(2).name; % e.g., 'M2'
    %    common_params.switch_parasitics.(switch_name_2) = struct('Cs', 2e-9, 'Rs', 0.2);
    % end
    % ... and so on for all N switches.

    % For buck_netlist_2.txt (assuming M1, M2 are identified in that order):
    common_params.switch_parasitics.M1 = struct('Cs', 2e-9, 'Rs', 0.2);
    common_params.switch_parasitics.M2 = struct('Cs', 2e-9, 'Rs', 0.2);
    % common_params.switch_parasitics.M3 = struct('Cs', 2e-9, 'Rs', 0.2);
    % common_params.switch_parasitics.M4 = struct('Cs', 2e-9, 'Rs', 0.2);
    % common_params.switch_parasitics.M5 = struct('Cs', 2e-8, 'Rs', 0.2);
    % common_params.switch_parasitics.M6 = struct('Cs', 2e-8, 'Rs', 0.2);
    % common_params.switch_parasitics.M7 = struct('Cs', 2e-9, 'Rs', 0.2);
    % common_params.switch_parasitics.M8 = struct('Cs', 2e-9, 'Rs', 0.2);
    fprintf('  (Example parasitics for M1, M2 are set. Modify for your netlist.)\n');
else
    fprintf('INFO: Parasitics are NOT included (common_params.include_parasitics is false).\n');
end

% --- Define Timing/Performance Analysis Options ---
common_params.num_timing_runs = 1; % Number of times to run the core calculation for timing stats
common_params.suppress_console_output = false; % Suppress detailed console output during runs

% --- Define Other Optional Solver/Optimization Parameters ---
% These are generally fine with defaults but can be tuned.
common_params.max_iterations = 50; % Max iterations for ZVS optimization
common_params.voltage_tolerance = 0.1; % Target voltage tolerance for ZVS (V)
common_params.damping_factor = 1;    % Damping factor for Newton-Raphson step in ZVS
common_params.min_dt = 1e-12;          % Absolute minimum deadtime allowed (s)
common_params.max_dt_fraction = 0.3;   % Old parameter, largely superseded by max_individual_dt_fraction_of_T
common_params.abstol = 1e-10;          % Absolute tolerance for calculateLambda Taylor series
common_params.degeneracy_threshold = 1e-20; % rcond threshold for detecting matrix degeneracy
common_params.significant_rate = 0;  % dV/dt threshold for ZVS heuristic (V/s)
common_params.adjust_increase = 1.01;  % Multiplicative factor for ZVS heuristic
common_params.adjust_decrease = 0.99;  % Multiplicative factor for ZVS heuristic
common_params.max_individual_dt_fraction_of_T = 0.25; % Max allowed individual DT as fraction of T_period
common_params.max_dt_step_change_fraction = 0.1; % Max change in DT per ZVS iteration, as fraction of T_period

% --- Tolerances for event_driven_timing_method (generateSequenceFromTimings) ---
common_params.time_tolerance_uniquetol = common_params.T_period * 1e-9; % For uniquetol in timeline generation
common_params.min_interval_duration = common_params.T_period * 1e-7;   % Min duration for an interval to be significant

% --- Configuration for Stage 1 (Parse Netlist, Get Valid States & Switch Info) ---
config_stage1 = struct();
config_stage1.Ron = common_params.Ron;
config_stage1.Roff = common_params.Roff;
config_stage1.fs = common_params.fs;
config_stage1.T_period = common_params.T_period;
config_stage1.include_parasitics = common_params.include_parasitics;
config_stage1.switch_parasitics = common_params.switch_parasitics; % Pass potentially populated struct
config_stage1.input_values = common_params.input_values;
config_stage1.stop_after_phase1 = true; % Crucial: Stop after identifying switches

% --- Run Stage 1 Analysis ---
% No need to edit anything within this try block below for Stage 1
fprintf('Calling MANSS for Phase 1 (Switch Identification)...\n');
try
    results_stage1 = MANSS(netlistFilePath, config_stage1);

    if isfield(results_stage1, 'switches') && isfield(results_stage1, 'valid_states')
        switches_identified_ordered = results_stage1.switches; % IMPORTANT: This defines the order
        valid_states_matrix_info = results_stage1.valid_states;
        num_valid_states_info = results_stage1.num_valid_states;

        fprintf('\n--- Stage 1 Complete --- \n');
        fprintf('Identified Switches (Order is CRITICAL for Section 2 configuration):\n');
        if isempty(switches_identified_ordered)
            fprintf('  No switches were identified in the netlist.\n');
        else
            for k_sw_info = 1:length(switches_identified_ordered)
                fprintf('  Switch %d: %s (Original Netlist Name)\n', k_sw_info, switches_identified_ordered(k_sw_info).name);
            end
        end
        fprintf('Number of switches found: %d\n', length(switches_identified_ordered));
        
        disp('Valid States Matrix (rows are states, cols correspond to identified switches):');
        if isempty(valid_states_matrix_info) && length(switches_identified_ordered) > 0
            disp('  (No valid states found after filtering - check netlist/config if switches exist)');
        elseif isempty(valid_states_matrix_info) && isempty(switches_identified_ordered)
            disp('  (No switches, so a single static state is assumed if circuit is valid)');
        else
            disp(valid_states_matrix_info);
        end
        fprintf('Found %d valid switching states (excluding all-ON by default).\n', num_valid_states_info);

        fprintf('\nACTION: Review "switches_identified_ordered" in your workspace.\n');
        fprintf('ACTION: If include_parasitics is true, ensure common_params.switch_parasitics in THIS SCRIPT is correctly defined for ALL identified switches using their printed names and order.\n');
        fprintf('ACTION: Then, configure Section 2 for these switches and run Section 2.\n');
    else
        fprintf(2,'\nStage 1 did not produce the expected "switches" or "valid_states" field in results_stage1.\n');
        if isfield(results_stage1, 'error')
            fprintf(2,'Error from MANSS: %s\n', results_stage1.error);
        end
        clear switches_identified_ordered valid_states_matrix_info num_valid_states_info; % Clear to prevent issues in Sec 2
    end
catch ME_stage1
    fprintf(2, '\n==================== STAGE 1 FAILED ====================\n');
    fprintf(2, 'Error Message: %s\n', ME_stage1.message);
    fprintf(2, 'Error occurred in:\n');
    for k_err = 1:length(ME_stage1.stack)
        fprintf(2, '  File: %s, Function: %s, Line: %d\n', ME_stage1.stack(k_err).file, ME_stage1.stack(k_err).name, ME_stage1.stack(k_err).line);
    end
    fprintf(2, '=========================================================\n');
    clear switches_identified_ordered valid_states_matrix_info num_valid_states_info; % Clear to prevent issues in Sec 2
end
fprintf('\n--- End of Section 1 ---\n');


%% =========================================================================
% SECTION 2: Configure Ideal Timings, Switches, Deadtimes, and Options for Full Analysis
% =========================================================================
fprintf('\n--- Running Section 2: Configure Full Analysis Options ---\n');

if ~exist('switches_identified_ordered', 'var')
    error('RunnerScript:SwitchInfoMissingS2', 'Variable "switches_identified_ordered" not found. Run Section 1 first.');
end
num_switches_found = length(switches_identified_ordered);
if num_switches_found ~= 2
    error('Expected to find 2 switches (M1, M2) for the Cuk converter, but found %d.', num_switches_found);
end

% --- MODIFICATION: Define Cuk Converter Specific Parameters ---
T_p = common_params.T_period; % 40us
user_D = 0.6;                 % 60% duty cycle
T_on = user_D * T_p;          % On-time for the primary switch (24us)

% --- MODIFICATION: Define Ideal Per-Switch ON/OFF Instants for Cuk Converter ---
% The two switches M1 and M2 operate complementarily.
% M1 (switch 1) is ON for the first part of the period (D*T).
% M2 (switch 2) is ON for the second part of the period ((1-D)*T).
user_ts_switches_ideal = cell(1, num_switches_found);

% Assuming M1 is the first switch found, and M2 is the second.
% This is based on the order in the netlist file.
fprintf('Configuring ideal timings for a 2-switch Cuk converter.\n');
% Switch 1 (M1): ON from 0 to T_on
user_ts_switches_ideal{1} = [0, T_on];
% Switch 2 (M2): ON from T_on to T_p
user_ts_switches_ideal{2} = [T_on, T_p];

fprintf('  Ideal timings for Cuk Converter (D=%.2f):\n', user_D);
fprintf('    Switch %s (idx 1): ON from %.2f us to %.2f us\n', switches_identified_ordered(1).name, user_ts_switches_ideal{1}(1)*1e6, user_ts_switches_ideal{1}(2)*1e6);
fprintf('    Switch %s (idx 2): ON from %.2f us to %.2f us\n', switches_identified_ordered(2).name, user_ts_switches_ideal{2}(1)*1e6, user_ts_switches_ideal{2}(2)*1e6);


user_defined_state_index_sequence = [];
user_operational_state_vectors = {};

% --- Define Deadtimes to be Applied ---
% These are the initial guesses for the optimization.
default_applied_deadtime_s = 1e-9; % 1 ns starting guess
user_deadtimes_applied_to_switches = ones(1, num_switches_found) * default_applied_deadtime_s;
fprintf('Initial guess for deadtime is %.1f ns for each switch.\n', default_applied_deadtime_s*1e9);

% --- ZVS Optimization ---
% Enable optimization and specify the target switch.
user_run_optimization = true;
user_zvs_switches_to_target = {'M1','M2'}; % Target M1 for Zero-Voltage Switching
fprintf('ZVS optimization is ON. Targeting switch: %s\n', user_zvs_switches_to_target{1});


% Define Complementary Switch Pairs for Cuk Converter ---
% The analyzer needs to know that M1 and M2 form a complementary pair.
% user_complementary_switch_pairs = { {'M1', 'M2'}, {'M3', 'M4'}, {'M5', 'M6'}};
user_complementary_switch_pairs = { {'M1', 'M2'} };
fprintf('User-defined complementary switch pairs: (%s, %s)\n', user_complementary_switch_pairs{1}{1}, user_complementary_switch_pairs{1}{2});

% --- Output Specification ---
% 'all_states' is the most comprehensive choice for initial analysis.
% The states will be i_L1, i_L2, v_C1, v_C2 (or similar).
user_outputs = {'all_states'};

fprintf('User configuration for full analysis defined in workspace.\n');
fprintf('Run Section 3 to execute the full analysis.\n');
fprintf('\n--- End of Section 2 ---\n');


%% =========================================================================
% SECTION 3: Run Full Analysis
% =========================================================================
fprintf('\n--- Running Section 3: Execute Full Analysis ---\n');
required_vars_s3 = {'netlistFilePath', 'common_params', 'switches_identified_ordered', ...
                    'user_ts_switches_ideal', 'user_deadtimes_applied_to_switches', ...
                    'user_run_optimization', 'user_zvs_switches_to_target', 'user_outputs'};
missing_vars_s3 = {};
for i_var_check_s3 = 1:length(required_vars_s3)
    if ~exist(required_vars_s3{i_var_check_s3}, 'var')
        missing_vars_s3{end+1} = required_vars_s3{i_var_check_s3};
    end
end
if common_params.include_parasitics && (~isfield(common_params, 'switch_parasitics') || isempty(fieldnames(common_params.switch_parasitics)) && num_switches_found > 0)
    warning('RunnerScript:PotentialMissingParasiticsS3', ...
            'common_params.include_parasitics is true and switches were found, but common_params.switch_parasitics appears empty. Ensure it was defined correctly in Section 1 for all identified switches.');
end
if ~isempty(missing_vars_s3)
    error('RunnerScript:MissingVariablesS3', ...
          'Cannot run Section 3. Missing variables: %s. Please run Sections 1 and 2 first, ensuring variables are created in the workspace.', ...
          strjoin(missing_vars_s3, ', '));
end
config_full = common_params;
config_full.ts_switches = user_ts_switches_ideal;
config_full.deadtimes_applied = user_deadtimes_applied_to_switches;
config_full.run_optimization = user_run_optimization;
config_full.zvs_switches_to_target = user_zvs_switches_to_target;
config_full.outputs = user_outputs;
config_full.switches_ordered_list = switches_identified_ordered;
config_full.user_complementary_switch_pairs = { {'M1', 'M2'} }; 
config_full.user_defined_state_index_sequence = user_defined_state_index_sequence;
config_full.user_operational_state_vectors = user_operational_state_vectors;
config_full.state_sequence = [];
config_full.state_duration_proportions = containers.Map();
config_full.initial_deadtimes = user_deadtimes_applied_to_switches;
config_full.stop_after_phase1 = false;
fprintf('\nCalling MANSS for full analysis...\n');
try
    results_full = MANSS(netlistFilePath, config_full);
    fprintf('\n--- Full Analysis Complete --- \n');
    if isfield(results_full, 'error') && ~isempty(results_full.error)
        fprintf(2, 'Analysis completed with an error flag from MANSS: %s\n', results_full.error);
    else
        fprintf('Access results in the "results_full" struct.\n');
    end
    
    if isfield(results_full, 'timing_stats') && isfield(results_full.timing_stats, 'num_runs') && results_full.timing_stats.num_runs > 1
        fprintf('\n--- Timing Analysis Results (Core Calculation) ---\n');
        stats = results_full.timing_stats;
        fprintf('  Number of Runs: %d\n', stats.num_runs);
        fprintf('  Mean Time:      %.4f ms\n', stats.mean_time_ms);
        fprintf('  Std Dev Time:   %.4f ms\n', stats.std_dev_time_ms);
        fprintf('  Min Time:       %.4f ms\n', stats.min_time_ms);
        fprintf('  Max Time:       %.4f ms\n', stats.max_time_ms);
        fprintf('--------------------------------------------------\n');
    elseif common_params.num_timing_runs > 1
        fprintf('\n--- Timing Analysis Warning: Timing stats missing or num_runs <=1 in results_full struct despite common_params.num_timing_runs > 1 ---\n');
    end
catch ME_stage3
    fprintf(2, '\n==================== FULL ANALYSIS FAILED (Section 3) ====================\n');
    fprintf(2, 'Error Message: %s\n', ME_stage3.message);
    fprintf(2, '---------------------------------------------------------\n');
    fprintf(2, 'Error occurred in:\n');
    for k_err = 1:length(ME_stage3.stack)
        fprintf(2, '  File: %s, Function: %s, Line: %d\n', ME_stage3.stack(k_err).file, ME_stage3.stack(k_err).name, ME_stage3.stack(k_err).line);
    end
    fprintf(2, '===========================================================================\n');
end
fprintf('\n--- Analysis Runner Script Finished ---\n');