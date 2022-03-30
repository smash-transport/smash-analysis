project(afterburner NONE)

set(afterburner_sims "")
set(afterburner_analysis "")
set(afterburner_plots "")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Plots)
function(afterburner_test
         test_name   # name of folder
    )
    message(STATUS "Configuring afterburner/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})
    file(MAKE_DIRECTORY ${test_folder}/Merged_spectra)


    set(smash_afterburner_config "${CMAKE_CURRENT_SOURCE_DIR}/afterburner.yaml")

    execute_process (
    COMMAND bash -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/${test_name} && gunzip -o sampling0.gz"
    OUTPUT_VARIABLE outVar
    )
    set(num_runs "20")
    foreach(i RANGE 1 ${num_runs})
      file(MAKE_DIRECTORY ${test_folder}/${i})
      file(MAKE_DIRECTORY ${test_folder}/${i}/data)
      file(MAKE_DIRECTORY ${test_folder}/${i}/analysis)
      set(final_particle_list "${test_folder}/${i}/data/particles_binary.bin")
      list(APPEND final_particle_lists ${final_particle_list})
    # Run afterburner evolution
      add_custom_command(OUTPUT "${final_particle_list}"
      COMMAND "${CMAKE_BINARY_DIR}/smash"
          "-i" "${smash_afterburner_config}"
          "-c" "Modi: { List: { File_Directory: ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}} }"
          "-o" "${test_folder}/${i}/data"
          "-f"
          ">" "/dev/null"
      COMMENT "Running SMASH afterburner for ${test_name} (${i}/${num_runs})."
      )


    #----------------------------------------------------------------------------#
    # Analyse data
    #----------------------------------------------------------------------------#
    # Define filenames for particle spectra and create list of those
      set(spectra_fnames "${test_folder}/${i}/analysis/yspectra.txt" ;
                         "${test_folder}/${i}/analysis/mtspectra.txt" ;
                         "${test_folder}/${i}/analysis/ptspectra.txt" ;
                         "${test_folder}/${i}/analysis/v2spectra.txt" ;
                         "${test_folder}/${i}/analysis/meanmt0_midrapidity.txt" ;
                         "${test_folder}/${i}/analysis/meanpt_midrapidity.txt" ;
                         "${test_folder}/${i}/analysis/midrapidity_yield.txt" ;
                         "${test_folder}/${i}/analysis/total_multiplicity.txt")
      set(analysis_output "")
      foreach(j ${spectra_fnames})
          list(APPEND analysis_output "${j}")
      endforeach(j)

      list(APPEND analysis_outputs ${analysis_output})


      # Perform analysis
      add_custom_command(OUTPUT ${analysis_output}
              COMMAND "python3" "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra.py"
                      "--output_files" ${analysis_output}
                      "--input_files" "${final_particle_list}"
              WORKING_DIRECTORY "${test_folder}"
              COMMENT "Analyzing spectra for ${test_name}  (${i}/${num_runs})."
              DEPENDS "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra.py"
          )
    endforeach(i)

    set(analysis_merged_outputs "")
    set(merged_spectra_fnames "${test_folder}/Merged_spectra/yspectra.txt" ;
                       "${test_folder}/Merged_spectra/mtspectra.txt" ;
                       "${test_folder}/Merged_spectra/ptspectra.txt" ;
                       "${test_folder}/Merged_spectra/v2spectra.txt" ;
                       "${test_folder}/Merged_spectra/meanmt0_midrapidity.txt" ;
                       "${test_folder}/Merged_spectra/meanpt_midrapidity.txt" ;
                       "${test_folder}/Merged_spectra/midrapidity_yield.txt" ;
                       "${test_folder}/Merged_spectra/total_multiplicity.txt")
    foreach(j ${merged_spectra_fnames})
        list(APPEND analysis_merged_outputs ${j})
    endforeach(j)
    add_custom_command(OUTPUT ${analysis_merged_outputs}
        COMMAND "python3" "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra.py"
                "--merge"
                "--input_files" ${analysis_outputs}
                "--output_files" ${analysis_merged_outputs}
        DEPENDS ${analysis_outputs}
        WORKING_DIRECTORY "${test_folder}"
        COMMENT "Merge spectra for ${test_name}."
    )

    set(target afterburner_${test_name})
    add_custom_target(${target} ALL DEPENDS ${final_particle_lists} ${analysis_merged_outputs} ${plot})
    #sim target for both processes
    set(target afterburner_sims_${test_name})
    add_custom_target(${target} ALL DEPENDS ${final_particle_lists})

    # analysis target for both processes
    set(target afterburner_analysis_${test_name})
    add_custom_target(${target} ALL DEPENDS ${analysis_merged_outputs})

    set(all_data_files ${all_data_files};${final_particle_lists} PARENT_SCOPE)
    set(all_analysis_outputs ${all_analysis_outputs};${analysis_merged_outputs} PARENT_SCOPE)
endfunction()

set(test_folder ${CMAKE_CURRENT_BINARY_DIR})
set(plot_names "${test_folder}/Plots/mtspectra_afterburner_RHIC_111.pdf" ;
               "${test_folder}/Plots/mtspectra_afterburner_LHC_111.pdf" ;)


afterburner_test("RHIC") #200GeV
afterburner_test("LHC")#5TeV

set(plot "")
foreach(j ${plot_names})
   list(APPEND plot "${j}")
endforeach(j)
if (NOT WITH_EXP_DATA) # plot without experimental data
  add_custom_command(OUTPUT ${plot}
      COMMAND "python3" "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra_plotter.py"
              "--theory" ${all_analysis_outputs}
              "--config" "${CMAKE_CURRENT_BINARY_DIR}/RHIC/data/config.yaml"
              "--output_folder"  "${test_folder}/Plots"
              "--comp_prev_version"
      DEPENDS "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra_plotter.py" ${all_analysis_outputs}
      COMMENT "Plotting results for ${test_name}."
  )
else() # plot with experimental data
  add_custom_command(OUTPUT ${plot}
      COMMAND "python3"  "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra_plotter.py"
              "--theory" ${all_analysis_outputs}
              "--experiment" "${EXP_DATA}/multiplicity/*/*.exp"
                             "${EXP_DATA}/midy_yield/*/*.exp"
                             "${EXP_DATA}/mean_pT/*/*.exp"
                             "${EXP_DATA}/mean_mT/*/*.exp"
              "--experiment_spectra" "${EXP_DATA}/mt_spectrum/*/*.exp"
                                     "${EXP_DATA}/y_spectrum/*/*.exp"
              "--config" "${CMAKE_CURRENT_BINARY_DIR}/RHIC/1/data/config.yaml"
              "--comp_prev_version"
              "--output_folder"  "${test_folder}/Plots"
      DEPENDS "${CMAKE_SOURCE_DIR}/test/energy_scan/mult_and_spectra_plotter.py" ${all_analysis_outputs}
      COMMENT "Plotting results for ${test_name} with experimental data."
  )
endif(NOT WITH_EXP_DATA)


# Target for all simulations
add_custom_target(afterburner_sims
    DEPENDS ${all_data_files}
)

# Target for the analysis
add_custom_target(afterburner_analysis
    DEPENDS ${all_analysis_outputs}
)

# Target for plotting
add_custom_target(afterburner_plots DEPENDS ${plot})


add_custom_target(afterburner DEPENDS afterburner_sims afterburner_analysis afterburner_plots)
