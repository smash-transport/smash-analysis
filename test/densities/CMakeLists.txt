project(densities NONE)

set(dens_box_sims "")
set(dens_box_plots "")

function(run_one_test
         test_name
         yaml_op_raw
         var_list
    )
    message(STATUS "Configuring densities/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})

    # Convert string to a list.
    separate_arguments(var_list)

    # run SMASH for both collision criteria
    set(i 0)
    set(data_files "")
    foreach(X IN LISTS var_list)
        set(data_density_box "${test_folder}/data/${i}")
        execute_process(COMMAND "echo" "${yaml_op_raw}"
                        COMMAND "sed" "s/_X_/${X}/"
                        OUTPUT_VARIABLE yaml_opt
                        OUTPUT_STRIP_TRAILING_WHITESPACE)

        # run box
        set(data_file "${data_density_box}/collisions_binary.bin")
        list(APPEND data_files ${data_file})
        add_custom_command(OUTPUT "${data_file}" "${data_density_box}/config.yaml"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                      "-o" "${data_density_box}"
                      "-i" "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                      "-p" "${CMAKE_CURRENT_SOURCE_DIR}/particles_${test_name}.txt"
                      "-d" "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes_${test_name}.txt"
                      "-f"
                      "-q"
            DEPENDS "${CMAKE_BINARY_DIR}/smash"
                    "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                    "${CMAKE_CURRENT_SOURCE_DIR}/particles_${test_name}.txt"
                    "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes_${test_name}.txt"
            WORKING_DIRECTORY "${test_folder}"
            COMMENT "Running SMASH for densities/${test_name}=${X}."
        )

        # FIXME: redirect SMASH output
        math(EXPR i "${i} + 1")
    endforeach()

    set(all_data_files ${all_data_files};${data_files} PARENT_SCOPE)
    set(el_box_sims ${el_box_sims};${target} PARENT_SCOPE)

    # analysis results


    # plot


    # targets for each setup
    set(target densities_${test_name}_sims)
    add_custom_target(${target} ALL DEPENDS ${data_files})

    set(target densities_${test_name}_analysis)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_files}"
            "${CMAKE_CURRENT_SOURCE_DIR}/elastic_scatrate.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
    )

    set(target densities_${test_name}_plots)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_plot}"
            "${CMAKE_CURRENT_SOURCE_DIR}/plotting/Ncoll_plotting.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
    )

endfunction()

# Pi-Rho System
run_one_test("pi_rho" " " "0 1 2 3 4 5 6 7 8 9 10")

# Pi K N L System
run_one_test("pi_k_n_lambda" " " "0 1 2 3 4 5 6 7 8 9 10")

# target for all simulations
add_custom_target(density_box_sims
    DEPENDS ${all_data_files}
)

# target for analysis
add_custom_target(density_box_analysis
    DEPENDS ${all_output_files}
)

# target for plotting
add_custom_target(density_box_plots
    DEPENDS ${el_box_plots}
)
