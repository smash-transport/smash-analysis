project(FOPI_pions NONE)

include(FloatMath)

set(FOPI_energies "")
set(all_data_files "")

function(run_one_energy
         energy
         )
  message(STATUS "Configuring FOPI_pions/${energy}.")
  set(test_folder "${CMAKE_CURRENT_BINARY_DIR}/${energy}")
  file(MAKE_DIRECTORY ${test_folder})
  set(data_files "")
  set(num_runs "10")
  foreach(i RANGE 1 ${num_runs})
      set(data_folder "${test_folder}/data/${i}")
      set(data_file "${data_folder}/particles_binary.bin")
      list(APPEND data_files ${data_file})
      # run smash
      add_custom_command(OUTPUT "${data_file}" "${data_folder}/config.yaml"
          COMMAND "${CMAKE_BINARY_DIR}/smash"
              "-c" "Modi: {Collider: {E_Kin: ${energy}} }"
              "-o" "${data_folder}"
              "-i" "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
              "-f"
              "-q"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/config.yaml"
              "${CMAKE_BINARY_DIR}/smash"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Running SMASH for FOPI_pions/${energy} (${i}/${num_runs})."
      )
  endforeach()
  # target for each energy
  set(target FOPI_pions_${energy})
  add_custom_target(${target} ALL
      DEPENDS "${data_files}"
  )
  set(all_data_files ${all_data_files};${data_files} PARENT_SCOPE)
  set(FOPI_energies ${FOPI_energies};${target} PARENT_SCOPE)
endfunction()

foreach(n RANGE 0 11)                  # loop over energy
	float_math(E "round(0.4 + ${n}*0.1, 1)")  # calc. energy
  run_one_energy(${E})
endforeach()

# analysis
add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/pion_yields.txt"
    COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/yields.py" ">pion_yields.txt"
    DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/yields.py"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Analyzing data for FOPI_pions."
)

# plot
if (NOT WITH_EXP_DATA) # plot without experimental data
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/FOPI_pions.pdf"
      COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
              "--comp_prev_version"
      DEPENDS
          "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
          "${CMAKE_CURRENT_BINARY_DIR}/pion_yields.txt"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Plotting data for FOPI_pions."
  )
else() # plot with experimental data
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/FOPI_pions.pdf"
      COMMAND "python3" "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
      "--FOPI_multiplicity" "${EXP_DATA}/multiplicity/AuAu/FOPI_pions.exp"
      "--FOPI_ratio" "${EXP_DATA}/piplus_piminus_ratio/AuAu/FOPI.exp"
      "--comp_prev_version"
      DEPENDS
          "${CMAKE_CURRENT_SOURCE_DIR}/plot.py"
          "${CMAKE_CURRENT_BINARY_DIR}/pion_yields.txt"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Plotting data for FOPI_pions."
  )
endif(NOT WITH_EXP_DATA)

add_custom_target(FOPI_pions_sims
    DEPENDS ${all_data_files}
)

add_custom_target(FOPI_pions_plots
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/FOPI_pions.pdf"
)
