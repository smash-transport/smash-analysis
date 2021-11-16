project(elastic_box NONE)

set(el_box_sims "")
set(el_box_plots "")

function(run_one_test
         test_name
         yaml_op_raw
         var_list
    )
    message(STATUS "Configuring elastic_box/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})

    # Convert string to a list.
    separate_arguments(var_list)

    # run SMASH for both collision criteria
    set(i 0)
    set(multipl_opt "Modi: {Box: {Init_Multiplicities: {111: 200}}}")
    # The cell length is explicitly reduced for the stochastic criterion below
    # to match the geometric criterion case, where the maximum cross section is
    # automatically set as the elastic cross section for elastic boxes.
    set(stoch_crit_opt "Collision_Term: {Collision_Criterion: Stochastic}, Collision_Term: {Fixed_Min_Cell_Length: 0.6}")
    set(geom_crit_opt "Collision_Term: {Collision_Criterion: Geometric}")
    set(data_files "")
    foreach(X IN LISTS var_list)
        set(data_geom_folder "${test_folder}/data-geometric/${i}")
        set(data_stoch_folder "${test_folder}/data-stochastic/${i}")
        set(data_cov_folder "${test_folder}/data-covariant/${i}")
        execute_process(COMMAND "echo" "${yaml_op_raw}"
                        COMMAND "sed" "s/_X_/${X}/"
                        OUTPUT_VARIABLE yaml_opt
                        OUTPUT_STRIP_TRAILING_WHITESPACE)

        # run with covariant criterion
        set(data_file "${data_cov_folder}/collisions_binary.bin")
        list(APPEND data_files ${data_file})
        add_custom_command(OUTPUT "${data_file}" "${data_cov_folder}/config.yaml"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                      "-o" "${data_cov_folder}"
                      "-i" "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                      "-p" "${CMAKE_CURRENT_SOURCE_DIR}/particles_only_pi0.txt"
                      "-d" "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes_all_off.txt"
                      "-c" "\\{${multipl_opt}, ${yaml_opt}\\}"
                      "-f"
                      "-q"
            DEPENDS "${CMAKE_BINARY_DIR}/smash"
                    "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                    "${CMAKE_CURRENT_SOURCE_DIR}/particles_only_pi0.txt"
                    "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes_all_off.txt"
            WORKING_DIRECTORY "${test_folder}"
            COMMENT "Running SMASH with covariant criterion for elastic_box/${test_name}=${X}."
        )

        # run with stochastic criterion
        set(data_file "${data_stoch_folder}/collisions_binary.bin")
        list(APPEND data_files ${data_file})
        add_custom_command(OUTPUT "${data_file}" "${data_stoch_folder}/config.yaml"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                      "-o" "${data_stoch_folder}"
                      "-i" "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                      "-p" "${CMAKE_CURRENT_SOURCE_DIR}/particles_only_pi0.txt"
                      "-d" "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes_all_off.txt"
                      "-c" "\\{${multipl_opt}, ${stoch_crit_opt} ,${yaml_opt}\\}"
                      "-f"
                      "-q"
            DEPENDS "${CMAKE_BINARY_DIR}/smash"
                    "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                    "${CMAKE_CURRENT_SOURCE_DIR}/particles_only_pi0.txt"
                    "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes_all_off.txt"
            WORKING_DIRECTORY "${test_folder}"
            COMMENT "Running SMASH with stochastic criterion for elastic_box/${test_name}=${X}."
        )

        # run with geometric criterion
        set(data_file "${data_geom_folder}/collisions_binary.bin")
        list(APPEND data_files ${data_file})
        add_custom_command(OUTPUT "${data_file}" "${data_geom_folder}/config.yaml"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                      "-o" "${data_geom_folder}"
                      "-i" "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                      "-p" "${CMAKE_CURRENT_SOURCE_DIR}/particles_only_pi0.txt"
                      "-d" "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes_all_off.txt"
                      "-c" "\\{${multipl_opt}, ${geom_crit_opt} ,${yaml_opt}\\}"
                      "-f"
                      ">" "/dev/null"
            DEPENDS "${CMAKE_BINARY_DIR}/smash"
                    "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                    "${CMAKE_CURRENT_SOURCE_DIR}/particles_only_pi0.txt"
                    "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes_all_off.txt"
            WORKING_DIRECTORY "${test_folder}"
            COMMENT "Running SMASH with geometric criterion for elastic_box/${test_name}=${X}."
        )

        # FIXME: redirect SMASH output
        math(EXPR i "${i} + 1")
    endforeach()

    set(all_data_files ${all_data_files};${data_files} PARENT_SCOPE)
    set(el_box_sims ${el_box_sims};${target} PARENT_SCOPE)

    # analysis results for geometric criterion
    set(output_file_geometric "${test_folder}/scatrate-geometric.txt")
    list(APPEND output_files ${output_file_geometric})
    add_custom_command(OUTPUT "${test_folder}/scatrate-geometric.txt"
        COMMAND "python2"
                  "${CMAKE_CURRENT_SOURCE_DIR}/elastic_scatrate.py"
                  "data-geometric" "scatrate-geometric.txt"
        WORKING_DIRECTORY "${test_folder}"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/elastic_scatrate.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
        COMMENT "Calculating scattering rate for geometric criterion for elastic_box/${test_name}."
    )

    # analysis results for stochastic criterion
    set(output_file_stochastic "${test_folder}/scatrate-stochastic.txt")
    list(APPEND output_files ${output_file_stochastic})
    add_custom_command(OUTPUT "${test_folder}/scatrate-stochastic.txt"
        COMMAND "python2"
                  "${CMAKE_CURRENT_SOURCE_DIR}/elastic_scatrate.py"
                  "data-stochastic" "scatrate-stochastic.txt"
        WORKING_DIRECTORY "${test_folder}"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/elastic_scatrate.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
        COMMENT "Calculating scattering rate for stochastic criterion for elastic_box/${test_name}."
    )

    # analysis results for covariant criterion
    set(output_file_covariant "${test_folder}/scatrate-covariant.txt")
    list(APPEND output_files ${output_file_covariant})
    add_custom_command(OUTPUT "${test_folder}/scatrate-covariant.txt"
        COMMAND "python2"
                  "${CMAKE_CURRENT_SOURCE_DIR}/elastic_scatrate.py"
                  "data-covariant" "scatrate-covariant.txt"
        WORKING_DIRECTORY "${test_folder}"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/elastic_scatrate.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
        COMMENT "Calculating scattering rate for covariant criterion for elastic_box/${test_name}."
    )

    set(all_output_files ${all_output_files};${output_files} PARENT_SCOPE)

    # plot
    set(output_plot "${test_folder}/scatrate.pdf")
    set(el_box_plots "${el_box_plots};${output_plot}" PARENT_SCOPE)
    add_custom_command(
        OUTPUT "${test_folder}/scatrate.pdf"
               "${test_folder}/${test_name}.txt"
        COMMAND "python2"
                  "${CMAKE_CURRENT_SOURCE_DIR}/plotting/Ncoll_plotting.py"
                  "--SMASH_data" "'${test_folder}/scatrate-*.txt'"
                  "--output" "${test_folder}/scatrate.pdf"
                  "--setup" "${test_name}"
                  "--comp_prev_version"
        WORKING_DIRECTORY "${test_folder}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/plotting/Ncoll_plotting.py"
                "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
        COMMENT "Plotting scattering rate for elastic_box/${test_name}."
    )

    # targets for each setup
    set(target elastic_box_${test_name}_sims)
    add_custom_target(${target} ALL DEPENDS ${data_files})

    set(target elastic_box_${test_name}_analysis)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_files}"
            "${CMAKE_CURRENT_SOURCE_DIR}/elastic_scatrate.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
    )

    set(target elastic_box_${test_name}_plots)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_plot}"
            "${CMAKE_CURRENT_SOURCE_DIR}/plotting/Ncoll_plotting.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
    )

endfunction()

# SCATTERING RATE VERSUS BOX VOLUME
run_one_test("scatrate_vs_V" "Modi: {Box: {Length: _X_}}" "3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 20.0")
#^ We might have to escape the braces by inserting "\\" before them.
#  It seems to work if we only escape the outmost braces in `run_one_test`.

# SCATTERING RATE VERSUS PARTICLE NUMBER
run_one_test("scatrate_vs_N" "Modi: {Box: {Init_Multiplicities: {111: _X_}}}" "20 40 60 100 120 140 160 180 200 240 280 330 400")

# SCATTERING RATE VERSUS CROSS-SECTION
run_one_test("scatrate_vs_sigma" "Collision_Term: {Elastic_Cross_Section: _X_}" "1.0 3.0 5.0 8.0 12.0 20.0 30.0 50.0 80.0")

# SCATTERING RATE VERSUS TEMPERATURE
run_one_test("scatrate_vs_T" "Modi: {Box: {Temperature: _X_}}" "0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9")

# SCATTERING RATE VERSUS TESTPARTICLES
run_one_test("scatrate_vs_Ntest" "General: {Testparticles: _X_}" "1 2 3 4 5 6 7 8 9 10")

# SCATTERING RATE VERSUS TIMESTEP
run_one_test("scatrate_vs_dt" "General: {Delta_Time: _X_}" "1.0 0.8 0.5 0.2 0.1 0.07 0.03 0.01 0.008 0.006 0.004")

# target for all simulations
add_custom_target(elastic_box_sims
    DEPENDS ${all_data_files}
)

# target for analysis
add_custom_target(elastic_box_analysis
    DEPENDS ${all_output_files}
)

# target for plotting
add_custom_target(elastic_box_plots
    DEPENDS ${el_box_plots}
)
