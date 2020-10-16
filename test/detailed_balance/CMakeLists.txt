project(detailed_balance NONE)

set(det_bal_sims "")
set(det_bal_analysis "")
set(det_bal_plots "")

function(one_balance_test_general
         test_name
         pdg_list
         reactions_to_analyze
         tstart_reactions_analysis
         reaction_counting_code
    )
    message(STATUS "Configuring detailed_balance/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})

    # run several instances of SMASH in parallel (with the same input)
    set(particles_files "")
    set(reactions_files "")
    foreach(i RANGE 1 20)
        set(data_folder "${test_folder}/data/${i}")
        set(collisions_file "${data_folder}/collisions_binary.bin")
        set(particles_file "${data_folder}/particles_binary.bin")
        list(APPEND particles_files ${particles_file})
        set(config_file "${data_folder}/config.yaml")
        add_custom_command(
            OUTPUT
                "${collisions_file}" "${particles_file}" "${config_file}"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                "-o" "${data_folder}"
                "-i" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/config.yaml"
                "-p" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/particles.txt"
                "-d" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/decaymodes.txt"
                "-f"
                ">" "/dev/null"
            DEPENDS "${CMAKE_BINARY_DIR}/smash"
                "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/config.yaml"
                "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/particles.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/decaymodes.txt"
            WORKING_DIRECTORY "${test_folder}"
            COMMENT "Running SMASH for detailed_balance/${test_name} (${i}/20)."
            VERBATIM
        )

        set(all_data_files ${all_data_files};${particles_files} PARENT_SCOPE)
        set(det_bal_sims ${det_bal_sims};${target} PARENT_SCOPE)

        # analyse number of reactions for each file
        set(data_folder "${test_folder}/data/${i}")
        set(reactions_file "${data_folder}/reaction_num_${i}.txt")
        list(APPEND reactions_files ${reactions_file})
        add_custom_command(OUTPUT ${reactions_file}
            COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/${reaction_counting_code}"
                "${reactions_file}" "\"${reactions_to_analyze}\""
                "${tstart_reactions_analysis}" "${config_file}" "${data_folder}/collisions_binary.bin"
            DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/${reaction_counting_code}"
                "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
                "${data_folder}/collisions_binary.bin"
            WORKING_DIRECTORY "${test_folder}"
            COMMENT "Analyzing number of reactions for detailed_balance/${test_name} (${i}/20)."
        )
    endforeach()

    # analyse multiplicities versus time
    set(multiplicity_file "${test_folder}/mult.txt")
    add_custom_command(OUTPUT "${test_folder}/mult.txt"
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/multiplicity_vs_time.py"
            "${test_folder}/mult.txt" "\"${pdg_list}\"" "${config_file}" ${particles_files}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/multiplicity_vs_time.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
        WORKING_DIRECTORY "${test_folder}"
        COMMENT "Calculating multiplicities for detailed_balance/${test_name}."
    )

    # merge output with reactions numbers to one file
    set(reaction_num_file "${test_folder}/reaction_num.txt")
    add_custom_command(OUTPUT "${test_folder}/reaction_num.txt"
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plotting/merge_reactionnum_output.py"
            "${test_folder}/reaction_num.txt" ${reactions_files}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/plotting/merge_reactionnum_output.py"
            "${reactions_files}"
        WORKING_DIRECTORY "${test_folder}"
        COMMENT "Merging numbers of reactions for detailed_balance/${test_name} to one file."
    )

    set(det_bal_analysis "${det_bal_analysis};${reactions_files};${multiplicity_file};${reaction_num_file}" PARENT_SCOPE)

    # plot
    set(output_plot_mult "${test_folder}/mult.pdf")
    add_custom_command(OUTPUT "${test_folder}/mult.pdf"
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plotting/plot_multiplicities_vs_time.py"
            "${test_folder}/mult.txt" "${test_folder}/mult.pdf" "${tstart_reactions_analysis}" "${config_file}"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/plotting/plot_multiplicities_vs_time.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
        WORKING_DIRECTORY "${test_folder}"
        COMMENT "Plotting multiplicities for detailed_balance/${test_name}."
    )
    set(output_plot_reacnum "${test_folder}/reaction_num.pdf")
    add_custom_command(OUTPUT "${test_folder}/reaction_num.pdf"
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plotting/plot_reactions_number.py"
            "--SMASH_data" "${test_folder}/reaction_num.txt"
            "--output" "${test_folder}/reaction_num.pdf"
            "--config" "${config_file}"
            "--setup" "${test_name}"
            "--comp_prev_version"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/plotting/plot_reactions_number.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
        WORKING_DIRECTORY "${test_folder}"
        COMMENT "Plotting reaction numbers for detailed_balance/${test_name}."
    )

    set(det_bal_plots "${det_bal_plots};${output_plot_mult};${output_plot_reacnum}" PARENT_SCOPE)

    # targets for each setup
    set(target detailed_balance_${test_name}_sims)
    add_custom_target(${target} ALL DEPENDS ${particles_files})

    set(target detailed_balance_${test_name}_analysis)
    add_custom_target(${target} ALL
        DEPENDS
            "${reactions_files}"
            "${multiplicity_file}"
            "${reaction_num_file}"
            "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
            "${CMAKE_CURRENT_SOURCE_DIR}/multiplicity_vs_time.py"
            "${CMAKE_CURRENT_SOURCE_DIR}/plotting/merge_reactionnum_output.py"
    )

    set(target detailed_balance_${test_name}_plots)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_plot_mult}"
            "${output_plot_reacnum}"
            "${CMAKE_CURRENT_SOURCE_DIR}/plotting/plot_multiplicities_vs_time.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
            "${CMAKE_CURRENT_SOURCE_DIR}/plotting/plot_reactions_number.py"
    )

endfunction()

function(one_balance_test
         test_name
         pdg_list
         reactions_to_analyze
         tstart_reactions_analysis
    )
  one_balance_test_general(${test_name} ${pdg_list} ${reactions_to_analyze} ${tstart_reactions_analysis}
                          "count_reactions_binQ2.py")
  set(all_data_files ${all_data_files} PARENT_SCOPE)
  set(det_bal_analysis ${det_bal_analysis} PARENT_SCOPE)
  set(det_bal_plots ${det_bal_plots} PARENT_SCOPE)
endfunction()

function(one_balance_test_angular
         test_name
         pdg_list
         reactions_to_analyze
         tstart_reactions_analysis
    )
  one_balance_test_general(${test_name} ${pdg_list} ${reactions_to_analyze} ${tstart_reactions_analysis}
                          "count_reactions_bin_t.py")
  set(all_data_files ${all_data_files} PARENT_SCOPE)
  set(det_bal_analysis ${det_bal_analysis} PARENT_SCOPE)
  set(det_bal_plots ${det_bal_plots} PARENT_SCOPE)
endfunction()

########## detailed balance for angular distributions, binning in -t Mandelstam variable #####

one_balance_test_angular("angular_NN_NDelta"
    "N⁺,N⁰,Δ⁺⁺,Δ⁺,Δ⁰,Δ⁻"
    "N⁺,N⁺:N⁰,Δ⁺⁺|3xN⁺,N⁺:N⁺,Δ⁺|3/2xN⁺,N⁰:N⁰,Δ⁺|3/2xN⁺,N⁰:N⁺,Δ⁰|3xN⁰,N⁰:N⁰,Δ⁰|N⁰,N⁰:N⁺,Δ⁻|N⁺,N⁺:Δ⁰,Δ⁺⁺|3/2xN⁺,N⁺:Δ⁺,Δ⁺|18/67xN⁺,N⁰:Δ⁻,Δ⁺⁺|18/43xN⁺,N⁰:Δ⁺,Δ⁰|3/2xN⁰,N⁰:Δ⁰,Δ⁰|N⁰,N⁰:Δ⁺,Δ⁻"
    "100.0"
)

########## thermal box, all possible reactions on, show reaction violating detailed balance most #####

one_balance_test("thermal_box"
    "π⁺,K⁺,N⁺,Λ,Σ⁻,Ξ⁻,ω,ρ⁺,Δ⁺,σ"
    "n_most_violating"
    "1800.0"
)

########## strangeness exchange ########

one_balance_test("Strangeness"
    "π⁺,π⁻,K⁺,K̅⁻,K⁰,K̅⁰,N⁺,N⁰,Δ⁺⁺,Δ⁺,Δ⁰,Δ⁻,Λ,Ξ⁻,Ω⁻"
    "n_most_violating"
    "2000.0"
)

one_balance_test("KN_KDelta"
    "K⁺,K̅⁻,K⁰,K̅⁰,N⁺,N⁰,Δ⁺⁺,Δ⁺,Δ⁰,Δ⁻"
    "3xK⁺,N⁺:K⁺,Δ⁺|K⁺,N⁺:K⁰,Δ⁺⁺|3xK⁺,N⁰:K⁺,Δ⁰|3xK⁺,N⁰:K⁰,Δ⁺|3xK⁰,N⁺:K⁰,Δ⁺|3xK⁰,N⁺:K⁺,Δ⁰|3xK⁰,N⁰:K⁰,Δ⁰|K⁰,N⁰:K⁺,Δ⁻|3xK̅⁻,N⁰:K̅⁻,Δ⁰|K̅⁻,N⁰:K̅⁰,Δ⁻|3xK̅⁻,N⁺:K̅⁰,Δ⁰|3xK̅⁻,N⁺:K̅⁻,Δ⁺"
    "40.0"
)


one_balance_test("N_pi_deutron"
    "π⁺,π⁻,π⁰,N⁺,N⁰,d,d',N̅⁰,N̅⁻,d̅,d̅'"
    "N⁺,N⁰:d'|N⁺,d:N⁺,d'|N⁰,d:N⁰,d'|π⁺,d:N⁺,N⁺|π⁰,d:N⁺,N⁰|π⁻,d:N⁰,N⁰|π⁺,d:π⁺,d'|π⁰,d:π⁰,d'|π⁻,d:π⁻,d'|N̅⁻,N̅⁰:d̅'|N̅⁻,d:N̅⁻,d'|N̅⁰,d:N̅⁰,d'|π⁺,d̅:N̅⁰,N̅⁰|π⁰,d̅:N̅⁻,N̅⁰|d̅,N⁰:d̅',N⁰|d̅,N̅⁰:d̅',N̅⁰"
    "500.0"
)

########## mesonic boxes ##########

one_balance_test("pi_rho"
    "π⁰,π⁺,π⁻,ρ⁰,ρ⁺,ρ⁻"
    "π⁰,π⁺:ρ⁺|π⁻,π⁰:ρ⁻|π⁻,π⁺:ρ⁰"
    "20.0"
)

one_balance_test("pi_rho_omega"
    "π⁰,π⁺,π⁻,ρ⁰,ρ⁺,ρ⁻,ω"
    "π⁰,ρ⁰:ω|π⁺,ρ⁻:ω|π⁻,ρ⁺:ω"
    "150.0"
)

one_balance_test("pi_sigma"
    "π⁰,π⁺,π⁻,σ"
    "2xπ⁰,π⁰:σ|π⁻,π⁺:σ"
    "20.0"
)

one_balance_test("pi_rho_sigma"
    "π⁰,π⁺,π⁻,ρ⁰,ρ⁺,ρ⁻,σ"
    "π⁰,π⁺:ρ⁺|π⁻,π⁰:ρ⁻|π⁻,π⁺:ρ⁰|2xπ⁰,π⁰:σ|π⁻,π⁺:σ"
    "20.0"
)

one_balance_test("pi_rho_f2"
    "π⁰,π⁺,π⁻,ρ⁰,ρ⁺,ρ⁻,f₂"
    "π⁰,π⁺:ρ⁺|π⁻,π⁰:ρ⁻|π⁻,π⁺:ρ⁰|π⁺,π⁻:f₂|2xπ⁰,π⁰:f₂|ρ⁺,ρ⁻:f₂|2xρ⁰,ρ⁰:f₂"
    "300.0"
)

########## multi-particle reactions ########

one_balance_test("multi_pi3_omega"
    "π⁰,π⁺,π⁻,ω"
    "π⁰,π⁺,π⁻:ω"
    "120.0"
)

one_balance_test("multi_pi3_phi"
    "π⁰,π⁺,π⁻,φ"
    "π⁰,π⁺,π⁻:φ"
    "120.0"
 )

one_balance_test("multi_pi2_eta_etaprime"
    "π⁰,π⁺,π⁻,η,η'"
    "π⁺,π⁻,η:η'|π⁰,π⁰,η:η'"
    "2500.0"
 )


########## baryonic boxes ##########

# N-pi-Delta, 1 <--> 2 reactions only (Delta <--> pi + N)
one_balance_test("N_pi_Delta_12"
    "π⁰,π⁺,π⁻,N⁺,N⁰,Δ⁺⁺,Δ⁺,Δ⁰,Δ⁻"
    "N⁺,π⁺:Δ⁺⁺|3/2xN⁺,π⁰:Δ⁺|3xN⁰,π⁺:Δ⁺|3xN⁺,π⁻:Δ⁰|3/2xN⁰,π⁰:Δ⁰|N⁰,π⁻:Δ⁻"
    "10.0"
)

# N-pi-Delta, 2 <--> 2 reactions only (N + N <--> N + Delta)
one_balance_test("N_pi_Delta_22"
    "π⁰,π⁺,π⁻,N⁺,N⁰,Δ⁺⁺,Δ⁺,Δ⁰,Δ⁻"
    "N⁺,N⁺:N⁰,Δ⁺⁺|3xN⁺,N⁺:N⁺,Δ⁺|3/2xN⁺,N⁰:N⁰,Δ⁺|3/2xN⁺,N⁰:N⁺,Δ⁰|3xN⁰,N⁰:N⁰,Δ⁰|N⁰,N⁰:N⁺,Δ⁻|N⁺,N⁺:Δ⁰,Δ⁺⁺|3/2xN⁺,N⁺:Δ⁺,Δ⁺|18/67xN⁺,N⁰:Δ⁻,Δ⁺⁺|18/43xN⁺,N⁰:Δ⁺,Δ⁰|3/2xN⁰,N⁰:Δ⁰,Δ⁰|N⁰,N⁰:Δ⁺,Δ⁻"
    "50.0"
)

# N-pi-Delta, all reactions (1<-->2 and 2<-->2)
one_balance_test("N_pi_Delta_all"
    "π⁰,π⁺,π⁻,N⁺,N⁰,Δ⁺⁺,Δ⁺,Δ⁰,Δ⁻"
    "N⁺,π⁺:Δ⁺⁺|3/2xN⁺,π⁰:Δ⁺|3xN⁰,π⁺:Δ⁺|3xN⁺,π⁻:Δ⁰|3/2xN⁰,π⁰:Δ⁰|N⁰,π⁻:Δ⁻|N⁺,N⁺:N⁰,Δ⁺⁺|3xN⁺,N⁺:N⁺,Δ⁺|3/2xN⁺,N⁰:N⁰,Δ⁺|3/2xN⁺,N⁰:N⁺,Δ⁰|3xN⁰,N⁰:N⁰,Δ⁰|N⁰,N⁰:N⁺,Δ⁻|N⁺,N⁺:Δ⁰,Δ⁺⁺|3/2xN⁺,N⁺:Δ⁺,Δ⁺|18/67xN⁺,N⁰:Δ⁻,Δ⁺⁺|18/43xN⁺,N⁰:Δ⁺,Δ⁰|3/2xN⁰,N⁰:Δ⁰,Δ⁰|N⁰,N⁰:Δ⁺,Δ⁻"
    "50.0"
)

# N(1440)
one_balance_test("N_1440"
    "π⁰,π⁺,π⁻,σ,N⁺,N⁰,Δ⁺⁺,Δ⁺,Δ⁰,Δ⁻,N(1440)⁺,N(1440)⁰"
    "2xN(1440)⁺:N⁺,π⁰|N(1440)⁺:N⁰,π⁺|2xN(1440)⁰:N⁰,π⁰|N(1440)⁰:N⁺,π⁻|N(1440)⁺:N⁺,σ|N(1440)⁰:N⁰,σ|N(1440)⁺,N⁺:N⁺,N⁺|N(1440)⁺,N⁰:N⁺,N⁰|N(1440)⁰,N⁺:N⁰,N⁺|N(1440)⁰,N⁰:N⁰,N⁰"
    "400.0"
)

#one_balance_test("N_K_Lambda"
#    "K⁺,K̅⁻,K⁰,K̅⁰,N⁺,N⁰,N(1650)⁺,N(1650)⁰,N(1720)⁺,N(1720)⁰,Λ,Λ(1520),Λ(1690),Λ(1820),Λ(1830),Λ(1890)"
#    "N⁺,K̅⁻:Λ(1520)|N⁺,K̅⁻:Λ(1690)|N⁺,K̅⁻:Λ(1820)|N⁺,K̅⁻:Λ(1830)|N⁺,K̅⁻:Λ(1890)|Λ,K⁺:N(1650)⁺|Λ,K⁺:N(1720)⁺"
#    "200.0"
#)

# target for all simulations
add_custom_target(detailed_balance_sims
    DEPENDS ${all_data_files}
)

# target for analysis
add_custom_target(detailed_balance_analysis
    DEPENDS ${det_bal_analysis}
)

# target for plotting
add_custom_target(detailed_balance_plots
    DEPENDS ${det_bal_plots}
)
