project(densities NONE)

set(dens_box_sims "")
set(dens_collider_sims "")
set(dens_plots "")

function(run_box_test
         test_name
         NFolder
    )
    message(STATUS "Configuring densities/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})

    # run SMASH box
    set(data_files "")
    foreach(X RANGE 0 ${NFolder}) # inclusive, runs until $NFolder
        set(data_density_box "${test_folder}/data/${X}")

        # run box
        set(data_file "${data_density_box}/collisions_binary.bin")
        list(APPEND data_files ${data_file})
        add_custom_command(OUTPUT "${data_file}" "${data_density_box}/config.yaml"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                      "-o" "${data_density_box}"
                      "-i" "${CMAKE_CURRENT_SOURCE_DIR}/box_setup/config.yaml"
                      "-p" "${CMAKE_CURRENT_SOURCE_DIR}/box_setup/particles_${test_name}.txt"
                      "-d" "${CMAKE_CURRENT_SOURCE_DIR}/box_setup/decaymodes_${test_name}.txt"
                      "-f"
                      "-q"
            DEPENDS "${CMAKE_BINARY_DIR}/smash"
                    "${CMAKE_CURRENT_SOURCE_DIR}/box_setup/config.yaml"
                    "${CMAKE_CURRENT_SOURCE_DIR}/box_setup/particles_${test_name}.txt"
                    "${CMAKE_CURRENT_SOURCE_DIR}/box_setup/decaymodes_${test_name}.txt"
            WORKING_DIRECTORY "${test_folder}"
            COMMENT "Running SMASH for densities/${test_name} (${X}/${NFolder})"
        )
    endforeach()

    set(all_data_files ${all_data_files};${data_files} PARENT_SCOPE)
    set(dens_box_sims ${dens_box_sims};${target} PARENT_SCOPE)

    # analysis results
    set(output_file_density "${test_folder}/density_grid.dat")
    list(APPEND output_files ${output_file_density})
    add_custom_command(
        OUTPUT "${test_folder}/density_grid.dat"
        COMMAND "python3"
                  "${CMAKE_CURRENT_SOURCE_DIR}/CalcDensityEvolution.py"
                  "${test_folder}/data/"
                  "${CMAKE_CURRENT_SOURCE_DIR}/box_setup/config.yaml"
                  "${test_folder}/density_grid.dat"

        WORKING_DIRECTORY "${test_folder}"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/CalcDensityEvolution.py"
        COMMENT "Calculating mean density grid for ${test_name}."
    )

    set(all_output_files ${all_output_files};${output_files} PARENT_SCOPE)

    # plot
    set(output_plot "${test_folder}/2d_density_${test_name}.pdf")
    set(dens_plots "${dens_plots};${output_plot}" PARENT_SCOPE)
    add_custom_command(
        OUTPUT "${test_folder}/2d_density_${test_name}.pdf"
        COMMAND "python3"
                  "${CMAKE_CURRENT_SOURCE_DIR}/plot_2d_density.py"
                  "${test_folder}/density_grid.dat"
                  "${test_folder}/2d_density_${test_name}.pdf"
        WORKING_DIRECTORY "${test_folder}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/plot_2d_density.py"
        COMMENT "Plotting 2d density densities/${test_name}."
    )

    # targets for each setup
    # sims
    set(target densities_${test_name}_sims)
    add_custom_target(${target} ALL DEPENDS ${data_files})

    # analysis
    set(target densities_${test_name}_analysis)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_files}"
            "${CMAKE_CURRENT_SOURCE_DIR}/CalcDensityEvolution.py"
    )

    # plotting
    set(target densities_${test_name}_plots)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_plot}"
            "${CMAKE_CURRENT_SOURCE_DIR}/plot_2d_density.py"

    )

endfunction()

function(run_collider_test
         test_name
         NFolder
    )
    message(STATUS "Configuring densities/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})

    # run SMASH collider
    set(data_files "")
    foreach(X RANGE 0 ${NFolder})
        set(data_density_collider "${test_folder}/data/${X}")

        set(data_file "${data_density_collider}/collisions_binary.bin")
        list(APPEND data_files ${data_file})
        add_custom_command(OUTPUT "${data_file}" "${data_density_collider}/config.yaml"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                      "-o" "${data_density_collider}"
                      "-i" "${CMAKE_CURRENT_SOURCE_DIR}/collider_setup/config.yaml"
                      "-f"
                      "-q"
            DEPENDS "${CMAKE_BINARY_DIR}/smash"
                    "${CMAKE_CURRENT_SOURCE_DIR}/collider_setup/config.yaml"
            WORKING_DIRECTORY "${test_folder}"
            COMMENT "Running SMASH for densities/${test_name} (${X}/${NFolder})"
        )
    endforeach()

    set(all_data_files ${all_data_files};${data_files} PARENT_SCOPE)
    set(dens_collider_sims ${dens_collider_sims};${target} PARENT_SCOPE)

    # analysis
    set(output_file_density "${test_folder}/density_central_cell.dat")
    list(APPEND output_files ${output_file_density})
    add_custom_command(
        OUTPUT "${test_folder}/density_central_cell.dat"
        COMMAND "python3"
                  "${CMAKE_CURRENT_SOURCE_DIR}/WriteOutDensityCentralCell.py"
                  "${test_folder}/data"
                  "${CMAKE_CURRENT_SOURCE_DIR}/collider_setup/config.yaml"
                  "${test_folder}/density_central_cell.dat"

        WORKING_DIRECTORY "${test_folder}"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/WriteOutDensityCentralCell.py"
        COMMENT "Analyse ${test_name}."
    )

    set(all_output_files ${all_output_files};${output_files} PARENT_SCOPE)

    # plot
    set(output_plot "${test_folder}/density_central_cell.pdf")
    set(dens_plots "${dens_plots};${output_plot}" PARENT_SCOPE)
    add_custom_command(
        OUTPUT "${test_folder}/density_central_cell.pdf"
        COMMAND "python3"
                  "${CMAKE_CURRENT_SOURCE_DIR}/plot_density_central_cell.py"
                  "${test_folder}/density_central_cell.dat"
                  "${test_folder}/density_central_cell.pdf"
                  "comp_prev_version"
        WORKING_DIRECTORY "${test_folder}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/plot_density_central_cell.py"
        COMMENT "Plotting central cell density densities/${test_name}."
    )


    # targets for each setup
    # sims
    set(target densities_${test_name}_sims)
    add_custom_target(${target} ALL DEPENDS ${data_files})

    # analysis
    set(target densities_${test_name}_analysis)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_files}"
            "${CMAKE_CURRENT_SOURCE_DIR}/WriteOutDensityCentralCell.py"
    )

    # plotting
    set(target densities_${test_name}_plots)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_plot}"
            "${CMAKE_CURRENT_SOURCE_DIR}/plot_density_central_cell.py"

    )


    ####################

endfunction()


# Pi-Rho System
run_box_test("pi_rho" "10")

# Pi K N L System
run_box_test("pi_k_n_lambda" "10")

# Au Au E_kin = 1.23
run_collider_test("AuAu_central_cell" "10")

# target for all simulations
add_custom_target(densities_sims
    DEPENDS ${all_data_files}
)

# target for analysis
add_custom_target(densities_analysis
    DEPENDS ${all_output_files}
)

# target for plotting
add_custom_target(densities_plots
    DEPENDS ${dens_plots}
)
