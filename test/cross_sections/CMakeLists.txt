project(cross_sections NONE)

include(FloatMath)

set(all_data_files "")
set(all_analysis_files "")
set(all_plots "")

# Create cross section targets.
function(cross_section_test
         test_name   # name of folder
         pdg1        # id of first particle
         pdg2        # id of second particle
         sqrts_min   # minimum energy
         sqrts_max   # maximum energy
         num_bins    # number of energy bins
         has_elastic # has data for the elastic cross section
    )
    message(STATUS "Configuring cross_sections/${test_name}.")
    float_math("dsqrts" "(${sqrts_max} - ${sqrts_min}) / ${num_bins}")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})

    set(analysis_files_1 "")
    set(analysis_files_2 "")
    set(analysis_files_3 "")
    set(analysis_files_4 "")
    set(analysis_files_5 "")
    set(analysis_files_6 "")
    set(analysis_files_7 "")
    set(data_files "")
    foreach(i RANGE 0 ${num_bins})
        float_math(energy "${sqrts_min} + ${i}*${dsqrts}")
        set(data_folder "${test_folder}/data/${i}")
        set(data_file "${data_folder}/collisions_binary.bin")
        list(APPEND data_files ${data_file})
        set(config_file "${data_folder}/config.yaml")

        # run SMASH
        add_custom_command(OUTPUT "${data_file}" "${config_file}"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                "-c" "Modi: {Collider: {Sqrtsnn: ${energy}, Projectile: {Particles: {${pdg1}: 1}}, Target: {Particles: {${pdg2}: 1}} } }"
                "-o" "${data_folder}"
                "-i" "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                "-f"
                "-q"
            DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
                "${CMAKE_BINARY_DIR}/smash"
            WORKING_DIRECTORY ${test_folder}
            COMMENT "Running SMASH for cross_sections/${test_name} at ${energy} GeV (${i}/${num_bins})."
            VERBATIM
        )

        set(all_data_files ${all_data_files};${data_files} PARENT_SCOPE)

        # analysis
        set(analysis_folder "${test_folder}/${i}")
        file(MAKE_DIRECTORY "${analysis_folder}")
        set(analysis_file_1 "${analysis_folder}/xs_individual.dat")
        set(analysis_file_2 "${analysis_folder}/xs_grouped.dat")
        set(analysis_file_3 "${analysis_folder}/xs_final_individual.dat")
        set(analysis_file_4 "${analysis_folder}/xs_final_grouped.dat")
        set(analysis_file_5 "${analysis_folder}/xs_final_individual_res.dat")
        set(analysis_file_6 "${analysis_folder}/xs_final_grouped_res.dat")
        set(analysis_file_7 "${analysis_folder}/xs_process_types.dat")
        list(APPEND analysis_files_1 ${analysis_file_1})
        list(APPEND analysis_files_2 ${analysis_file_2})
        list(APPEND analysis_files_3 ${analysis_file_3})
        list(APPEND analysis_files_4 ${analysis_file_4})
        list(APPEND analysis_files_5 ${analysis_file_5})
        list(APPEND analysis_files_6 ${analysis_file_6})
        list(APPEND analysis_files_7 ${analysis_file_7})
        set(all_analysis_files "${all_analysis_files};${analysis_files_1};${analysis_files_2};${analysis_files_3};${analysis_files_4};${analysis_files_5};${analysis_files_6};${analysis_files_7}" PARENT_SCOPE)
        add_custom_command(OUTPUT "${analysis_file_1}" "${analysis_file_2}" "${analysis_file_3}" "${analysis_file_4}" "${analysis_file_5}" "${analysis_file_6}" "${analysis_file_7}"
            VERBATIM
            COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/geometrical_cross_section.py"
                "--verbose"
                "${pdg1}" "${pdg2}" "${config_file}"
                "-o1" "${analysis_file_1}"
                "-o2" "${analysis_file_2}"
                "-o3" "${analysis_file_3}"
                "-o4" "${analysis_file_4}"
                "-o5" "${analysis_file_5}"
                "-o6" "${analysis_file_6}"
                "-o7" "${analysis_file_7}"
                "--resonances" "φ,Δ,ρ,η,η',ω"
                "${data_file}"
            DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/geometrical_cross_section.py"
                "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
            WORKING_DIRECTORY ${test_folder}
            COMMENT "Analyzing data for cross_sections/${test_name} (${i}/${num_bins})."
        )
    endforeach()

    add_custom_command(OUTPUT "${test_folder}/xs_individual.dat"
        VERBATIM
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
            "${test_folder}/xs_individual.dat"
            "--data" ${analysis_files_1}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Combining results for Analyzing data for cross_sections/${test_name} (individual)."
    )
    add_custom_command(OUTPUT "${test_folder}/xs_grouped.dat"
        VERBATIM
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
            "${test_folder}/xs_grouped.dat"
            "--data" ${analysis_files_2}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Combining results for Analyzing data for cross_sections/${test_name} (grouped)."
    )
    add_custom_command(OUTPUT "${test_folder}/xs_final_individual.dat"
        VERBATIM
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
            "${test_folder}/xs_final_individual.dat"
            "--data" ${analysis_files_3}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Combining results for Analyzing data for cross_sections/${test_name} (final individual)."
    )
    add_custom_command(OUTPUT "${test_folder}/xs_final_grouped.dat"
        VERBATIM
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
            "${test_folder}/xs_final_grouped.dat"
            "--data" ${analysis_files_4}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Combining results for Analyzing data for cross_sections/${test_name} (final grouped)."
    )
    add_custom_command(OUTPUT "${test_folder}/xs_final_individual_res.dat"
        VERBATIM
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
            "${test_folder}/xs_final_individual_res.dat"
            "--data" ${analysis_files_5}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Combining results for Analyzing data for cross_sections/${test_name} (final individual for resonance)."
    )
    add_custom_command(OUTPUT "${test_folder}/xs_final_grouped_res.dat"
        VERBATIM
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
            "${test_folder}/xs_final_grouped_res.dat"
            "--data" ${analysis_files_6}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Combining results for Analyzing data for cross_sections/${test_name} (final grouped for resonance)."
    )
    add_custom_command(OUTPUT "${test_folder}/xs_process_type.dat"
        VERBATIM
        COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
            "${test_folder}/xs_process_type.dat"
            "--data" ${analysis_files_7}
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/combine.py"
            ${analysis_files_7}
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Combining results for Analyzing data for cross_sections/${test_name} (process type)."
    )

    # plotting
    if(has_elastic AND WITH_EXP_DATA)
        set(data_elastic "${EXP_DATA}/cross_sections/${test_name}/data_elastic.exp")
        set(data_elastic_arg "--elastic_xs")
    else()
        set(data_elastic "")
        set(data_elastic_arg "")
    endif()

    if (NOT WITH_EXP_DATA) # plot without experimental data
      set(output_plot_1 "${test_folder}/xs_individual.pdf")
      add_custom_command(OUTPUT ${test_folder}/xs_individual.pdf
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "xs_individual.dat" "${test_folder}/xs_individual.pdf"
              "${data_elastic_arg}" "${data_elastic}" "${test_name}"
              "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${test_folder}/xs_individual.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (individual) for cross_sections/${test_name}."
      )
      set(output_plot_2 "${test_folder}/xs_grouped.pdf")
      add_custom_command(OUTPUT ${test_folder}/xs_grouped.pdf
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "xs_grouped.dat" "${test_folder}/xs_grouped.pdf"
              "${data_elastic_arg}" "${data_elastic}" "${test_name}"
              "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${test_folder}/xs_grouped.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (grouped) for cross_sections/${test_name}."
      )
      set(output_plot_3 "${test_folder}/xs_final_grouped.pdf")
      add_custom_command(OUTPUT ${test_folder}/xs_final_grouped.pdf
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "xs_final_grouped.dat" "${test_folder}/xs_final_grouped.pdf"
              "${data_elastic_arg}" "${data_elastic}" "${test_name}"
              "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${test_folder}/xs_final_grouped.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (final grouped) for cross_sections/${test_name}."
      )
      set(output_plot_4 "${test_folder}/xs_process_type.pdf")
      add_custom_command(OUTPUT ${test_folder}/xs_process_type.pdf
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "xs_process_type.dat" "${test_folder}/xs_process_type.pdf"
              "${data_elastic_arg}" "${data_elastic}" "${test_name}"
              "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${test_folder}/xs_process_type.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (process type) for cross_sections/${test_name}."
      )
    else() # plot with experimental data
      set(output_plot_1 "${test_folder}/xs_individual.pdf")
      add_custom_command(OUTPUT ${test_folder}/xs_individual.pdf
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "xs_individual.dat" "${test_folder}/xs_individual.pdf"
              "--total_xs" "${EXP_DATA}/cross_sections/${test_name}/data_total.exp"
              "${data_elastic_arg}" "${data_elastic}" "${test_name}"
              "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${EXP_DATA}/cross_sections/${test_name}/data_total.exp"
              "${test_folder}/xs_individual.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (individual) for cross_sections/${test_name}."
      )
      set(output_plot_2 "${test_folder}/xs_grouped.pdf")
      add_custom_command(OUTPUT ${test_folder}/xs_grouped.pdf
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "xs_grouped.dat" "${test_folder}/xs_grouped.pdf"
              "--total_xs" "${EXP_DATA}/cross_sections/${test_name}/data_total.exp"
              "${data_elastic_arg}" "${data_elastic}" "${test_name}"
              "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${EXP_DATA}/cross_sections/${test_name}/data_total.exp"
              "${test_folder}/xs_grouped.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (grouped) for cross_sections/${test_name}."
      )
      set(output_plot_3 "${test_folder}/xs_final_grouped.pdf")
      add_custom_command(OUTPUT ${test_folder}/xs_final_grouped.pdf
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "xs_final_grouped.dat" "${test_folder}/xs_final_grouped.pdf"
              "--total_xs" "${EXP_DATA}/cross_sections/${test_name}/data_total.exp"
              "${data_elastic_arg}" "${data_elastic}" "${test_name}"
              "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${EXP_DATA}/cross_sections/${test_name}/data_total.exp"
              "${test_folder}/xs_final_grouped.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (final grouped) for cross_sections/${test_name}."
      )
      set(output_plot_4 "${test_folder}/xs_process_type.pdf")
      add_custom_command(OUTPUT ${test_folder}/xs_process_type.pdf
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "xs_process_type.dat" "${test_folder}/xs_process_type.pdf"
              "--total_xs" "${EXP_DATA}/cross_sections/${test_name}/data_total.exp"
              "${data_elastic_arg}" "${data_elastic}" "${test_name}"
              "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${EXP_DATA}/cross_sections/${test_name}/data_total.exp"
              "${test_folder}/xs_process_type.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (process type) for cross_sections/${test_name}."
      )
    endif(NOT WITH_EXP_DATA)

    set(all_plots "${all_plots};${output_plot_1};${output_plot_2};${output_plot_3};${output_plot_4}" PARENT_SCOPE)

    # targets for each cross section
    set(target cross_sections_${test_name}_sims)
    add_custom_target(${target} ALL DEPENDS ${data_files})

    set(target cross_sections_${test_name}_analysis)
    add_custom_target(${target} ALL
        DEPENDS
            "${analysis_files_1}"
            "${analysis_files_2}"
            "${analysis_files_3}"
            "${analysis_files_4}"
            "${analysis_files_5}"
            "${analysis_files_6}"
            "${analysis_files_7}"
            "${CMAKE_CURRENT_SOURCE_DIR}/geometrical_cross_section.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/smash_basic_scripts.py"
    )

    set(target cross_sections_${test_name}_plots)
    add_custom_target(${target} ALL
        DEPENDS
            "${output_plot_1}"
            "${output_plot_2}"
            "${output_plot_3}"
            "${output_plot_4}"
            ${cross_section_plots_${test_name}}
            "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
            "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
    )

endfunction()

# Plot partial cross sections.
#
# Has to be called before `cross_section_test`.
function(partial_cross_section
         test_name   # name of folder
         data_name   # name of experimental data
         final_state # only the given final state is considered
         is_final_state_unstable # whether final state is unstable
    )
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    set(output ${test_folder}/xs_${data_name}.pdf)
    if(is_final_state_unstable)
      set(analyzed_result "xs_final_individual_res.dat")
    else()
      set(analyzed_result "xs_final_individual.dat")
    endif()

    if (NOT WITH_EXP_DATA) # plot without experimental data
      add_custom_command(OUTPUT ${output}
          VERBATIM
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "--final_state" "${final_state}"
              "${analyzed_result}" "${output}" "${test_name}"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${test_folder}/xs_final_individual.dat"
              "${test_folder}/xs_final_individual_res.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (${data_name}) for cross_sections/${test_name}."
      )
    else() # plot with experimental data
      add_custom_command(OUTPUT ${output}
          VERBATIM
          COMMAND "python2" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "--final_state" "${final_state}"
              "--final_state_data" "${EXP_DATA}/cross_sections/${test_name}/${data_name}.exp"
              "${analyzed_result}" "${output}" "${test_name}"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "${CMAKE_SOURCE_DIR}/python_scripts/common_plotting.py"
              "${EXP_DATA}/cross_sections/${test_name}/${data_name}.exp"
              "${test_folder}/xs_final_individual.dat"
              "${test_folder}/xs_final_individual_res.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data (${data_name}) for cross_sections/${test_name}."
      )
    endif(NOT WITH_EXP_DATA)

    set(cross_section_plots_${test_name} ${cross_section_plots_${test_name}};${output} PARENT_SCOPE)
    set(all_plots ${all_plots};${output} PARENT_SCOPE)
endfunction()


# If you add a partial cross section for an unstable particle, make sure to add
# it to the resonances argument `--resonances` for the analysis! Otherwise the
# final state resonance will not be reconstructed by the analysis.

cross_section_test("piplus_piminus" "211" "-211" "0.28" "1.58" "65" "yes")

#partial_cross_section("piminus_proton" "lambda_pizero_kzero" "Λ+π⁰+K⁰" "no")  # not seen in SMASH
#partial_cross_section("piminus_proton" "lambda_piminus_kplus" "Λ+π⁻+K⁺" "no")   # not seen in SMASH
partial_cross_section("piminus_proton" "lambda_X" "Λ+X" "no")
#partial_cross_section("piminus_proton" "lambda_pizero_piminus_kplus" "Λ+π⁰+π⁻+K⁺" "no")  # not seen in SMASH
#partial_cross_section("piminus_proton" "lambda_piplus_piminus_kzero" "Λ+π⁺+π⁻+K⁰" "no")  # not seen in SMASH
partial_cross_section("piminus_proton" "lambda_kzero" "Λ+K⁰" "no")
partial_cross_section("piminus_proton" "sigmaminus_kplus" "Σ⁻+K⁺" "no")
cross_section_test("piminus_proton" "-211" "2212" "1.08" "5.08" "300" "yes")

partial_cross_section("piplus_proton" "sigmaplus_kplus" "Σ⁺+K⁺" "no")
cross_section_test("piplus_proton"  "211" "2212" "1.08" "5.08" "200" "yes")

cross_section_test("antiproton_proton" "-2212" "2212" "1.96" "6.96" "200" "yes") 

partial_cross_section("proton_proton" "neutron_deltaplusplus" "n+Δ⁺⁺" "yes")
partial_cross_section("proton_proton" "piplus_X" "π⁺+X" "no")
partial_cross_section("proton_proton" "proton_proton_pizero" "p+p+π⁰" "no")
partial_cross_section("proton_proton" "proton_neutron_piplus" "p+n+π⁺" "no")
partial_cross_section("proton_proton" "proton_proton_pizero_pizero" "p+p+π⁰+π⁰" "no")
partial_cross_section("proton_proton" "proton_proton_piplus_piminus" "p+p+π⁺+π⁻" "no")
partial_cross_section("proton_proton" "proton_neutron_piplus_pizero" "p+n+π⁺+π⁰" "no")
partial_cross_section("proton_proton" "lambda_X" "Λ+X" "no")
partial_cross_section("proton_proton" "lambda_proton_kplus" "Λ+p+K⁺" "no")
partial_cross_section("proton_proton" "kplus_X" "K⁺+X" "no")
partial_cross_section("proton_proton" "kzero_X" "K⁰+X" "no")
partial_cross_section("proton_proton" "rhozero_X" "ρ⁰+X" "yes")
partial_cross_section("proton_proton" "proton_proton_rhozero" "p+p+ρ⁰" "yes")
#partial_cross_section("proton_proton" "lambda_sigmazero" "Λ+Σ⁰" "no")  # charge not conserved in Landolt-Börnstein!
partial_cross_section("proton_proton" "proton_proton_kzero_kzero" "p+p+K⁰+K̅⁰" "no")  # extremely small cross section
partial_cross_section("proton_proton" "sigmazero_proton_kplus" "Σ⁰+p+K⁺" "no")
partial_cross_section("proton_proton" "sigmazero_X" "Σ⁰+X" "no")
partial_cross_section("proton_proton" "sigmaplus_neutron_kplus" "Σ⁺+n+K⁺" "no")
partial_cross_section("proton_proton" "sigmaplus_proton_kzero" "Σ⁺+p+K⁰" "no")
#partial_cross_section("proton_proton" "phi_X" "φ+X" "yes")  # extremely small cross section
partial_cross_section("proton_proton" "proton_proton_eta" "p+p+η" "yes")
#partial_cross_section("proton_proton" "proton_proton_etaprime" "p+p+η'" "yes")  # not seen in SMASH
partial_cross_section("proton_proton" "proton_proton_omegaMeson" "p+p+ω" "yes")
#partial_cross_section("proton_proton" "proton_proton_kplus_kminus" "p+p+K⁺+K̅⁻" "no")  # extremely small cross section
#partial_cross_section("proton_proton" "proton_proton_kplus_kminus_via_phi" "p+p+φ" "yes")  # extremely small cross section
cross_section_test("proton_proton" "2212" "2212" "1.96" "8.96" "140" "yes")

partial_cross_section("neutron_proton" "neutron_proton_eta" "n+p+η" "yes")
partial_cross_section("neutron_proton" "neutron_proton_pizero" "n+p+π⁰" "no")
partial_cross_section("neutron_proton" "proton_proton_piminus" "p+p+π⁻" "no")
partial_cross_section("neutron_proton" "neutron_proton_piplus_piminus" "n+p+π⁺+π⁻" "no")
partial_cross_section("neutron_proton" "proton_proton_pizero_piminus" "p+p+π⁰+π⁻" "no")
cross_section_test("neutron_proton" "2112" "2212" "1.96" "8.96" "140" "yes")
# cross_section_test("proton_neutron" "2212" "2112" "1.96" "4.96" "60" "yes")

partial_cross_section("kminus_proton" "lambda_pizero" "Λ+π⁰" "no")
partial_cross_section("kminus_proton" "sigmaminus_piplus" "Σ⁻+π⁺" "no")
partial_cross_section("kminus_proton" "sigmazero_pizero" "Σ⁰+π⁰" "no")
partial_cross_section("kminus_proton" "sigmaplus_piminus" "Σ⁺+π⁻" "no")
cross_section_test("kminus_proton" "-321" "2212" "1.433" "2.532" "55" "yes")

partial_cross_section("kplus_proton" "deltaplusplus_kzero" "Δ⁺⁺+K⁰" "yes")
partial_cross_section("kplus_proton" "deltaplus_kplus" "Δ⁺+K⁺" "yes")
cross_section_test("kplus_proton" "321" "2212" "1.433" "2.532" "55" "yes")

partial_cross_section("kminus_neutron" "lambda_piminus" "Λ+π⁻" "no")
partial_cross_section("kminus_neutron" "sigmazero_piminus" "Σ⁰+π⁻" "no")
partial_cross_section("kminus_neutron" "sigmaminus_pizero" "Σ⁻+π⁰" "no")
partial_cross_section("kminus_neutron" "deltaminus_kzero" "Δ⁻+K̅⁰" "yes")
partial_cross_section("kminus_neutron" "deltazero_kminus" "Δ⁰+K̅⁻" "yes")
cross_section_test("kminus_neutron" "-321" "2112" "1.433" "2.532" "55" "yes")

partial_cross_section("kplus_neutron" "proton_piminus_kplus" "p+π⁻+K⁺" "no")
partial_cross_section("kplus_neutron" "proton_kzero" "p+K⁰" "no")
cross_section_test("kplus_neutron" "321" "2112" "1.433" "2.532" "55" "no")

# Target for all simulations
add_custom_target(cross_sections_sims
    DEPENDS ${all_data_files}
)

# Target for analysis
add_custom_target(cross_sections_analysis
    DEPENDS ${all_analysis_files}
)

# Target for plotting
add_custom_target(cross_sections_plots
    DEPENDS ${all_plots}
)
