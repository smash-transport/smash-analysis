project(angular_distributions NONE)

set(ang_dists_sims "")
set(ang_dists_analysis "")
set(ang_dists_plots "")

function(angular_test
         test_name   # name of folder
         pdg1        # id of first particle
         pdg2        # id of second particle
         plab        # projectile momentum in lab frame [GeV], target at rest
    )
    message(STATUS "Configuring angular_distributions/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})

    # run SMASH
    set(data_files "")
    foreach(i RANGE 1 50)
        set(data_folder "${test_folder}/data/${i}")
        set(data_file "${data_folder}/collisions_binary.bin")
        list(APPEND data_files ${data_file})
        set(config_file "${data_folder}/config.yaml")
        add_custom_command(OUTPUT "${data_file}" "${config_file}"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                "-c" "Modi: {Collider: {P_Lab: ${plab}, Projectile: {Particles: {${pdg1}: 1}}, Target: {Particles: {${pdg2}: 1}} } }"
                "-o" "${data_folder}"
                "-i" "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                "-f"
                "-q"
            DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                "${CMAKE_BINARY_DIR}/smash"
            WORKING_DIRECTORY ${test_folder}
            COMMENT "Running SMASH for angular_distributions/${test_name} (${i}/50)."
            VERBATIM
        )
    endforeach()

    set(all_data_files ${all_data_files};${data_files} PARENT_SCOPE)
    set(ang_dists_sims ${ang_dists_sims};${target} PARENT_SCOPE)


    # analysis
    set(output_theta "${test_folder}/theta.dat")
    set(output_t "${test_folder}/t.dat")
    set(ang_dists_analysis "${ang_dists_analysis};${output_t};${output_theta}" PARENT_SCOPE)
    add_custom_command(OUTPUT "${test_folder}/theta.dat" "${test_folder}/t.dat"
        COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/angular.py"
            "${pdg1}" "${pdg2}" "${config_file}" ${data_files}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/angular.py"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Analyzing data for angular_distributions/${test_name}."
    )

    # plot theta
    set(output_plot_theta "${test_folder}/theta.pdf")
    add_custom_command(OUTPUT ${test_folder}/theta.pdf
        COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/plot_theta.py"
                "--SMASH_data" "theta.dat"
                "--pdg1" "${pdg1}"
                "--pdg2" "${pdg2}"
                "--plab" "${plab}"
                "--setup" "${test_name}"
                "--comp_prev_version"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/plot_theta.py"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Plotting theta for angular_distributions/${test_name}."
    )

    # plot t
    set(output_plot_t "${test_folder}/t.pdf")
    if (NOT WITH_EXP_DATA) # plot without experimental data
      add_custom_command(OUTPUT ${test_folder}/t.pdf
          COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/plot_t.py"
                  "--SMASH_data" "t.dat"
                  "--pdg1" "${pdg1}"
                  "--pdg2" "${pdg2}"
                  "--plab" "${plab}"
                  "--setup" "${test_name}"
                  "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot_t.py"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting t for angular_distributions/${test_name}."
      )
    else() # plot with experimental data
      add_custom_command(OUTPUT ${test_folder}/t.pdf
          COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/plot_t.py"
                  "--SMASH_data" "t.dat"
                  "--pdg1" "${pdg1}"
                  "--pdg2" "${pdg2}"
                  "--plab" "${plab}"
                  "--setup" "${test_name}"
                  "--comp_prev_version"
                  "--exp_data" "${EXP_DATA}/angular_distributions/"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot_t.py"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting t for angular_distributions/${test_name}."
      )
    endif(NOT WITH_EXP_DATA)

    set(ang_dists_plots "${ang_dists_plots};${output_plot_theta};${output_plot_t}" PARENT_SCOPE)

    # targets for each setup
    set(target angular_distributions_${test_name}_sims)
    add_custom_target(${target} ALL DEPENDS ${data_files})

    set(target angular_distributions_${test_name}_analysis)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_theta}"
            "${output_t}"
            "${CMAKE_CURRENT_SOURCE_DIR}/angular.py"
    )

    set(target angular_distributions_${test_name}_plots)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_plot_theta}"
            "${output_plot_t}"
            "${CMAKE_CURRENT_SOURCE_DIR}/plot_theta.py"
            "${CMAKE_CURRENT_SOURCE_DIR}/plot_t.py"
    )

endfunction()

# p+p
angular_test("pp_1.25" "2212" "2212" "1.249")
angular_test("pp_1.75" "2212" "2212" "1.75")
angular_test("pp_2.80" "2212" "2212" "2.8")

# n+p
angular_test("np_1.194" "2112" "2212" "1.194")

# Target for all simulations
add_custom_target(angular_distributions_sims
    DEPENDS ${all_data_files}
)

# Target for the analysis
add_custom_target(angular_distributions_analysis
    DEPENDS ${ang_dists_analysis}
)

# Target for plotting
add_custom_target(angular_distributions_plots
    DEPENDS ${ang_dists_plots}
)
