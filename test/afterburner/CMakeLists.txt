project(afterburner NONE)

set(afterburner_sims "")
set(afterburner_analysis "")
set(afterburner_plots "")


function(afterburner_test
         test_name   # name of folder
    )
    message(STATUS "Configuring afterburner/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})

    set(smash_afterburner_config "${CMAKE_CURRENT_SOURCE_DIR}/afterburner.yaml")
    set(final_particle_list "${test_folder}/particles_binary.bin")

    list(APPEND final_particle_lists ${final_particle_list})
    set(all_data_files ${all_data_files};${final_particle_lists} PARENT_SCOPE)

    execute_process (
    COMMAND bash -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/${test_name} && unzip -o sampling0.zip"
    OUTPUT_VARIABLE outVar
    )

    # Run afterburner evolution
    add_custom_command(OUTPUT "${final_particle_list}"
    COMMAND "${CMAKE_BINARY_DIR}/smash"
        "-i" "${smash_afterburner_config}"
        "-c" "Modi: { List: { File_Directory: ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}} }"
        "-o" "${test_folder}"
        "-f"
        ">" "/dev/null"
    COMMENT "Running SMASH afterburner for ${test_name}."
    )
    set(target afterburner_${test_name})

    #----------------------------------------------------------------------------#
    # Analyse data
    #----------------------------------------------------------------------------#
    # Define filenames for particle spectra and create list of those
    set(spectra_fnames "${test_folder}/yspectra.txt" ;
                       "${test_folder}/mtspectra.txt" ;
                       "${test_folder}/ptspectra.txt" ;
                       "${test_folder}/v2spectra.txt" ;
                       "${test_folder}/meanmt0_midrapidity.txt" ;
                       "${test_folder}/meanpt_midrapidity.txt" ;
                       "${test_folder}/midrapidity_yield.txt" ;
                       "${test_folder}/total_multiplicity.txt")
    set(analysis_outputs "")
    foreach(j ${spectra_fnames})
        list(APPEND analysis_output "${j}")
    endforeach(j)

    list(APPEND analysis_outputs ${analysis_output})
    set(all_analysis_outputs ${all_analysis_outputs};${analysis_outputs} PARENT_SCOPE)

    # Perform analysis
    add_custom_command(OUTPUT ${analysis_output}
            COMMAND "python2" "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra.py"
                    "--output_files" ${analysis_outputs}
                    "--input_files" "${final_particle_list}"
            COMMENT "Analyzing spectra for $for ${test_name})."
            DEPENDS "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra.py"
        )
    add_custom_target(${target} ALL DEPENDS ${final_particle_list} ${analysis_output})
    #----------------------------------------------------------------------------#
    # Plot spectra
    #----------------------------------------------------------------------------#
    # Define plot names
    set(plot_names "${test_folder}/yspectra.pdf" ;
                   "${test_folder}/mtspectra.pdf" ;
                   "${test_folder}/ptspectra.pdf" ;
                   "${test_folder}/v2spectra.pdf" ;
                   "${test_folder}/meanmt0_midrapidity.pdf" ;
                   "${test_folder}/meanpt_midrapidity.pdf" ;
                   "${test_folder}/midrapidity_yield.pdf" ;
                   "${test_folder}/total_multiplicity.pdf")

    set(plot "")
    foreach(j ${plot_names})
       list(APPEND plot "${j}")
    endforeach(j)

    set(mockup "${test_folder}/afterburner_mockup.pdf")

    # Perform plotting
    if (NOT WITH_EXP_DATA) # plot without experimental data
      add_custom_command(OUTPUT "${mockup}"
          COMMAND "python" "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra_plotter.py"
                  "--theory" ${all_analysis_outputs}
                  "--config" "${CMAKE_CURRENT_BINARY_DIR}/${test_name}/config.yaml"
                  "--comp_prev_version"
          DEPENDS ${all_analysis_outputs}
          COMMENT "Plotting results for ${target_name}."
      )
    #[[else() # plot with experimental data
      add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/afterburner_mockup.pdf"
          COMMAND "python" "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra_plotter.py"
                  "--theory" ${all_analysis_outputs}
                  "--experiment" "${EXP_DATA}/multiplicity/*/*.exp"
                                 "${EXP_DATA}/midy_yield/*/*.exp"
                                 "${EXP_DATA}/mean_pT/*/*.exp"
                                 "${EXP_DATA}/mean_mT/*/*.exp"
                  "--experiment_spectra" "${EXP_DATA}/mt_spectrum/*/*.exp"
                                         "${EXP_DATA}/y_spectrum/*/*.exp"
                  "--config" "${CMAKE_CURRENT_BINARY_DIR}/${test_name}/config.yaml"
                  "--comp_prev_version"
          DEPENDS ${all_analysis_outputs}
          COMMENT "Plotting results for ${target_name}."
      )#]]
    endif()
endfunction()

afterburner_test("RHIC") #200GeV
afterburner_test("LHC")#5TeV

# Target for all simulations
add_custom_target(afterburner_sims
    DEPENDS ${all_data_files}
)

# Target for the analysis
add_custom_target(afterburner_analysis
    DEPENDS ${all_analysis_outputs}
)

# Target for plotting
add_custom_target(afterburner_plots DEPENDS "${mockup}")


add_custom_target(afterburner DEPENDS afterburner_sims afterburner_analysis afterburner_plots)
