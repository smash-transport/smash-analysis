project(energy_scan NONE)

include(FloatMath)

set(all_data_files "")
set(all_analysis_outputs "")

function(run_one_energy
         energy_and_system
         )
  list(GET energy_and_system 0 energy)
  list(GET energy_and_system 1 system)
  message(STATUS "Configuring energy_scan/${energy}/${system}.")
  set(test_folder "${CMAKE_CURRENT_BINARY_DIR}/${energy}/${system}")
  file(MAKE_DIRECTORY ${test_folder})
  set(data_files "")
  set(analysis_outputs "")
  set(spectra_fnames "yspectra.txt" ; "mtspectra.txt" ; "ptspectra.txt" ; "v2spectra.txt" ;
                     "meanmt0_midrapidity.txt" ; "meanpt_midrapidity.txt" ;
                     "midrapidity_yield.txt" ; "total_multiplicity.txt")
  set(num_runs "20")
  foreach(i RANGE 1 ${num_runs})
      set(data_folder "${test_folder}/data/${i}")
      set(data_file "${data_folder}/particles_binary.bin")
      list(APPEND data_files ${data_file})
      # run smash
      add_custom_command(OUTPUT "${data_file}" "${data_folder}/config.yaml"
          COMMAND "${CMAKE_BINARY_DIR}/smash"
              "-c" "Modi: {Collider: {Sqrtsnn: ${energy}} }"
              "-o" "${data_folder}"
              "-i" "${CMAKE_CURRENT_SOURCE_DIR}/config_${system}.yaml"
              "-f"
              "-q"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/config_${system}.yaml"
              "${CMAKE_BINARY_DIR}/smash"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Running SMASH for energy_scan/${energy}/${system} (${i}/${num_runs})."
      )

      # analyze smash output, create files with spectra, multiplicities, mean pt, etc
      set(analysis_outputs_i "")
      foreach(j ${spectra_fnames})
          list(APPEND analysis_outputs_i "${data_folder}/${j}")
      endforeach(j)
      list(APPEND analysis_outputs ${analysis_outputs_i})
      add_custom_command(OUTPUT ${analysis_outputs_i}
          COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/mult_and_spectra.py"
                  "--output_files" ${analysis_outputs_i}
                  "--input_files" "${data_file}"
          WORKING_DIRECTORY "${test_folder}"
          COMMENT "Analyzing spectra for energy_scan/${energy}/${system} (${i}/${num_runs})."
      )
  endforeach(i)

  set(analysis_merged_outputs "")
  foreach(j ${spectra_fnames})
      list(APPEND analysis_merged_outputs "${test_folder}/${j}")
  endforeach(j)
  add_custom_command(OUTPUT ${analysis_merged_outputs}
      COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/mult_and_spectra.py"
              "--merge"
              "--input_files" ${analysis_outputs}
              "--output_files" ${analysis_merged_outputs}
      DEPENDS ${analysis_outputs}
      WORKING_DIRECTORY "${test_folder}"
      COMMENT "Merging analyzed spectra for energy_scan/${energy}/${system}."
  )

  # sims target for each energy
  set(target energy_scan_${energy}_${system}_sims)
  add_custom_target(${target} ALL DEPENDS ${data_files})

  # analysis target for each energy
  set(target energy_scan_${energy}_${system}_analysis)
  add_custom_target(${target} ALL DEPENDS ${analysis_merged_outputs})

  set(all_data_files ${all_data_files};${data_files} PARENT_SCOPE)
  set(all_analysis_outputs ${all_analysis_outputs};${analysis_merged_outputs} PARENT_SCOPE)
endfunction()


run_one_energy("2.695;pp")
run_one_energy("3.325;pp")
run_one_energy("3.845;pp")
run_one_energy("4.305;pp")
run_one_energy("4.915;pp")
run_one_energy("6.405;pp")
run_one_energy("7.735;pp")
run_one_energy("8.865;pp")
run_one_energy("12.395;pp")
run_one_energy("17.31;pp")
run_one_energy("32.0;pp")
run_one_energy("62.4;pp")
run_one_energy("130;pp")
run_one_energy("200;pp")
run_one_energy("900;pp")
run_one_energy("7000;pp")

run_one_energy("2.695;AuAu")
run_one_energy("3.325;AuAu")
run_one_energy("3.845;AuAu")
run_one_energy("4.305;AuAu")
run_one_energy("4.915;PbPb")
run_one_energy("6.405;PbPb")
run_one_energy("7.735;PbPb")
run_one_energy("8.865;PbPb")
run_one_energy("12.395;PbPb")
run_one_energy("17.31;PbPb")
run_one_energy("5.5;AuAu")
run_one_energy("7.7;AuAu")
run_one_energy("11.5;AuAu")
run_one_energy("19.6;AuAu")
run_one_energy("27;AuAu")
run_one_energy("39;AuAu")
run_one_energy("62.4;AuAu")
run_one_energy("130;AuAu")
run_one_energy("200;AuAu")
run_one_energy("2760;PbPb")
run_one_energy("5020;PbPb")



add_custom_target(energy_scan_sims DEPENDS ${all_data_files})
add_custom_target(energy_scan_analysis DEPENDS ${all_analysis_outputs})

add_custom_target(energy_scan_plots DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/mockup.pdf")

if (NOT WITH_EXP_DATA) # plot without experimental data
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/mockup.pdf"
      COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/mult_and_spectra_plotter.py"
              "--theory" ${all_analysis_outputs}
              "--config" "${CMAKE_CURRENT_BINARY_DIR}/2.695/pp/data/1/config.yaml"
              "--comp_prev_version"
      DEPENDS ${all_analysis_outputs}
      COMMENT "Plotting results for energy_scan without expdata. ${all_analysis_outputs}"
  )
else() # plot with experimental data
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/mockup.pdf"
      COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/mult_and_spectra_plotter.py"
              "--theory" ${all_analysis_outputs}
              "--experiment" "${EXP_DATA}/multiplicity/*/*.exp"
                             "${EXP_DATA}/midy_yield/*/*.exp"
                             "${EXP_DATA}/mean_pT/*/*.exp"
                             "${EXP_DATA}/mean_mT/*/*.exp"
              "--experiment_spectra" "${EXP_DATA}/mt_spectrum/*/*.exp"
                                     "${EXP_DATA}/y_spectrum/*/*.exp"
              "--config" "${CMAKE_CURRENT_BINARY_DIR}/2.695/pp/data/1/config.yaml"
              "--comp_prev_version"
      DEPENDS ${all_analysis_outputs}
      COMMENT "Plotting results for energy_scan with experimental data. ${EXP_DATA}"
  )
endif(NOT WITH_EXP_DATA)

add_custom_target(energy_scan DEPENDS energy_scan_sims energy_scan_analysis energy_scan_plots)
