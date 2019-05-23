project(dileptons NONE)

set(dilepton_tests "")

function(dilepton_test
         system    # collision system
         Ekin      # kinetic energy
         # optional args:
         # * acc_file -- file with acceptance matrix
         # * plep_min -- minimum lepton momentum
         # * plep_max -- maximum lepton momentum
        )

    if (ARGC GREATER 2 AND ARGV2)
      # do filtering!
      set(test_name "${system}_${Ekin}_filtered")
      set(acc_file "${EXP_DATA}/HADES_dilepton_filtering/${ARGV2}")
      set(plep_min "${ARGV3}")
      set(plep_max "${ARGV4}")
      set(data_dir "${EXP_DATA}")
      set(filter "${CMAKE_CURRENT_SOURCE_DIR}/HAFT.so")
      # pion average needed for CC and ArKCl
      if (system STREQUAL "CC" OR system STREQUAL "ArKCl")
        set(pion_average_total_file "${CMAKE_CURRENT_BINARY_DIR}/${test_name}/other.avg_pion.dat")
      else()
        set(pion_average_total_file "")
      endif()
      # set resolution
      if (system STREQUAL "CC")
        set(resolution "1")   # low res.
      else()
        set(resolution "3")   # high res.
      endif()
    else()
      # no filter
      set(test_name "${system}_${Ekin}_unfiltered")
      # no data
      set(data_dir "")
    endif()

    message(STATUS "Configuring dileptons/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    set(def_file "${CMAKE_CURRENT_SOURCE_DIR}/definitions.py")
    file(MAKE_DIRECTORY ${test_folder})

    # set config (CC needs specific imp. param. per energy)
    if (system STREQUAL "CC")
      set(config_file "${CMAKE_CURRENT_SOURCE_DIR}/${system}/config_${Ekin}.yaml")
      set(kinetic_energy_config "")
    else()
      set(config_file "${CMAKE_CURRENT_SOURCE_DIR}/${system}/config.yaml")
      set(kinetic_energy_config "-c Modi: {Collider: { E_Kin: ${Ekin} } }")
    endif()

    set(hist_files "")
    set(pion_files "")
    set(num_runs 100)
    foreach(i RANGE 1 ${num_runs})
        set(data_folder "${test_folder}/data/${i}")
        set(data_file "${data_folder}/Dileptons.bin")
        set(hist_file "${data_folder}/hist_mass.txt")
        set(pion_file "${data_folder}/other.avg_pion.dat")
        list(APPEND hist_files ${hist_file})
        list(APPEND pion_files ${pion_file})
        if (system STREQUAL "CC" OR system STREQUAL "ArKCl")
          set(output_files "${data_file}" "${data_folder}/particle_lists.oscar" ${data_folder}/config.yaml)
        else()
          set(output_files "${data_file}" ${data_folder}/config.yaml)
        endif()


        # run SMASH
        add_custom_command(OUTPUT ${output_files}
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                "-o" "${data_folder}"
                "-i" "${config_file}"
                "-d" "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes.txt"
                "-f"
                ${kinetic_energy_config}
                ">" "/dev/null"
            DEPENDS
                "${config_file}"
                "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes.txt"
                "${CMAKE_BINARY_DIR}/smash"
            WORKING_DIRECTORY ${test_folder}
            COMMENT "Running SMASH for dileptons/${test_name} at ${Ekin} GeV (${i}/${num_runs})."
            VERBATIM
            )

        # run analysis
        add_custom_command(OUTPUT "${hist_file}" "${data_folder}/other.version.dat"
            COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/create_hist_from_binary_output.py" ${data_file}
                    "${acc_file}" "${plep_min}" "${plep_max}" "${resolution}"
            DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/create_hist_from_binary_output.py"
                ${filter}
                ${data_file}
                ${def_file}
            WORKING_DIRECTORY ${data_folder}
            COMMENT "Analyzing data for dileptons/${test_name}  (${i}/${num_runs})."
            )

        # count also pions for CC or ArKCl plots
        add_custom_command(OUTPUT "${pion_file}"
            COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/avg_pions.py" "${data_folder}/particle_lists.oscar"
            DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/avg_pions.py"
                "${data_folder}/particle_lists.oscar"
            WORKING_DIRECTORY ${data_folder}
            COMMENT "Averaging pions for dileptons/${test_name}  (${i}/${num_runs})."
            )
    endforeach()

    # compile HAFT code with f2py
    add_custom_command(OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/HAFT.so"
        COMMAND "f2py" "-m" "HAFT" "-c" "${EXP_DATA}/HADES_dilepton_filtering/readHAFT2.f90"
                "skip:" "momspread" "readhaftmatrix" "getmatrixval" "gettableval"
                "param" "smearhades4momentum" "readhaftpairmatrix"
        DEPENDS "${EXP_DATA}/HADES_dilepton_filtering/readHAFT2.f90"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Running f2py."
    )

    # add up all histograms
    add_custom_command(OUTPUT "${test_folder}/hist_mass.txt" "${test_folder}/other.version.dat"
        COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/add_hist.py"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/add_hist.py"
            ${hist_files}
            ${def_file}
            "${test_folder}/data/1/other.version.dat"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Adding histograms for dileptons/${test_name}."
    )

    # adding pion averages for total average for CC or ArKCl
    add_custom_command(OUTPUT "${test_folder}/other.avg_pion.dat"
        COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/add_avg.py"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/add_avg.py"
            ${pion_files}
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Adding pion averages for dileptons/${test_name}."
    )

    # plotting
    if (NOT WITH_EXP_DATA) # plot without experimental data
      add_custom_command(OUTPUT ${test_folder}/plot_mass.pdf
          COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/plot_hist.py"
                  "--system" "${system}"
                  "--energy" "${Ekin}"
                  "--comp_prev_version"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot_hist.py"
              "${test_folder}/hist_mass.txt"
              ${pion_average_total_file}
              "${def_file}"
              "${test_folder}/other.version.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data for dileptons/${test_name}."
      )
    else() # plot with experimental data
      add_custom_command(OUTPUT ${test_folder}/plot_mass.pdf
          COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/plot_hist.py"
                  "--system" "${system}"
                  "--energy" "${Ekin}"
                  "--comp_prev_version"
                  "--data_dir" "${data_dir}"
          DEPENDS
              "${CMAKE_CURRENT_SOURCE_DIR}/plot_hist.py"
              "${test_folder}/hist_mass.txt"
              ${pion_average_total_file}
              "${def_file}"
              "${test_folder}/other.version.dat"
          WORKING_DIRECTORY ${test_folder}
          COMMENT "Plotting data for dileptons/${test_name}."
      )
    endif(NOT WITH_EXP_DATA)

    # add a target for each test case (i.e. particular system & energy)
    set(target dileptons_${test_name})
    add_custom_target(${target} ALL
                      DEPENDS "${test_folder}/plot_mass.pdf")
    # add a collective target for each system (including all energies)
    if (TARGET dileptons_${system})
      add_dependencies(dileptons_${system} ${target})
    else()
      add_custom_target(dileptons_${system}
                        DEPENDS ${target})
    endif()


    # add a target for analysis suite
    # tests to be included in analysis (all unfiltered):
    #   pp at 1.25 and 3.5 GeV
    #   CC at 2.0 GeV
    #   ArKCl at 1.76 GeV

    string(COMPARE EQUAL "${test_name}" pp_1.25_unfiltered pp_1.25)
    string(COMPARE EQUAL "${test_name}" pp_3.5_unfiltered pp_3.5)
    string(COMPARE EQUAL "${test_name}" CC_2.0_unfiltered CC2.0)
    string(COMPARE EQUAL "${test_name}" ArKCl_1.76_unfiltered ArKCl1.76)

    if (pp_1.25 OR pp_3.5 OR CC2.0 OR ArKCl1.76)
      set(target_ana dileptons_ana_${test_name})
      add_custom_target(${target_ana} ALL
                        DEPENDS "${test_folder}/plot_mass.pdf")
      add_custom_target(dileptons_ana_${system}_${Ekin}
                          DEPENDS ${target_ana})
    endif()

    set(dileptons_ana ${dileptons_ana};${target_ana} PARENT_SCOPE)

endfunction()


function(box_test
         system    # box system
        )

    set(test_name "box/${system}")
    message(STATUS "Configuring dileptons/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})
    set(def_file "${CMAKE_CURRENT_SOURCE_DIR}/definitions.py")


    set(hist_files "")
    set(num_runs 100)
    foreach(i RANGE 1 ${num_runs})
        set(data_folder "${test_folder}/data/${i}")
        set(data_file "${data_folder}/Dileptons.bin")
        set(hist_file "${data_folder}/hist_mass.txt")
        list(APPEND hist_files ${hist_file})
        # run SMASH
        add_custom_command(OUTPUT "${data_file}" "${data_folder}/config.yaml"
            COMMAND "${CMAKE_BINARY_DIR}/smash"
                "-o" "${data_folder}"
                "-i" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/config.yaml"
                "-p" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/particles.txt"
                "-d" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/decaymodes.txt"
                "-f"
            DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/config.yaml"
                "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/particles.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/decaymodes.txt"
                "${CMAKE_BINARY_DIR}/smash"
            WORKING_DIRECTORY ${test_folder}
            COMMENT "Running SMASH for dileptons/${test_name} (${i}/${num_runs})."
            VERBATIM
            )
        # run analysis
        add_custom_command(OUTPUT "${hist_file}"
            COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/create_hist_from_binary_output.py" ${data_file}
            DEPENDS
                "${CMAKE_CURRENT_SOURCE_DIR}/create_hist_from_binary_output.py"
                ${data_file}
                ${def_file}
            WORKING_DIRECTORY ${data_folder}
            COMMENT "Analyzing data for dileptons/${test_name}  (${i}/${num_runs})."
            )
    endforeach()

    # add up all histograms
    add_custom_command(OUTPUT "${test_folder}/hist_mass.txt"
        COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/add_hist.py"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/add_hist.py"
            ${hist_files}
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Adding histograms for dileptons/${test_name}."
    )

    # plotting
    add_custom_command(OUTPUT ${test_folder}/plot_mass.pdf
        COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/plot_hist.py"  "box" "0.2"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/plot_hist.py"
            "${test_folder}/hist_mass.txt"
            "${def_file}"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Plotting data for dileptons/${test_name}."
    )

    # add a target for each test case
    set(target dileptons_box_${system})
    add_custom_target(${target} ALL
                      DEPENDS "${test_folder}/plot_mass.pdf")
endfunction()


function(shining_BR_test
        )

    set(test_name "shining_BR_test")
    message(STATUS "Configuring dileptons/${test_name}.")
    set(test_folder ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    file(MAKE_DIRECTORY ${test_folder})

    set(data_folder "${test_folder}/data")
    set(data_dilepton_file "${data_folder}/Dileptons.oscar")
    set(data_hadrons_file "${data_folder}/full_event_history.oscar")
    set(hist_dilepton_file "${data_folder}/hist_mass.txt")
    set(hist_pion_file "${data_folder}/hist_mass_pi.txt")

    # run SMASH
    add_custom_command(OUTPUT "${data_dilepton_file}" "${data_hadrons_file}"
        COMMAND "${CMAKE_BINARY_DIR}/smash"
            "-o" "${data_folder}"
            "-d" "${CMAKE_CURRENT_SOURCE_DIR}/decaymodes.txt"
            "-i" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/config.yaml"
            "-f"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/config.yaml"
            "${CMAKE_BINARY_DIR}/smash"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Running SMASH for dileptons/${test_name}."
        VERBATIM
        )

    # run analysis dileptons
    add_custom_command(OUTPUT "${hist_dilepton_file}"
        COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/analysis.py" "${data_dilepton_file}" "11"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/analysis.py"
            "${data_dilepton_file}"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Analyzing dilepton output for dileptons/${test_name}."
        )

    # run analysis pions
    add_custom_command(OUTPUT "${hist_pion_file}"
        COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/analysis.py" "${data_hadrons_file}" "211"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/analysis.py"
            "${data_hadrons_file}"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Analyzing hadron output for dileptons/${test_name}."
        )

    # plotting multiplicities
    add_custom_command(OUTPUT "${test_folder}/plot_mass_n.pdf"
        COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/plot_multip.py"  "box" "0.2"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/plot_multip.py"
            "${hist_dilepton_file}"
            "${hist_pion_file}"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Plotting multiplicities for dileptons/${test_name}."
    )

    # plotting BR
    add_custom_command(OUTPUT "${test_folder}/plot_mass_br.pdf"
        COMMAND "python" "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/plot_br.py"  "box" "0.2"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/plot_br.py"
            "${hist_pion_file}"
            "${hist_dilepton_file}"
        WORKING_DIRECTORY ${test_folder}
        COMMENT "Plotting BR for dileptons/${test_name}."
    )

    # add a target for each test case
    set(target dileptons_shining_BR_test)
    add_custom_target(${target} ALL
                      DEPENDS "${test_folder}/plot_mass_n.pdf" "${test_folder}/plot_mass_br.pdf")
endfunction()




## tests without filtering ##

dilepton_test("pp" "1.25")
dilepton_test("pp" "2.2")
dilepton_test("pp" "3.5")
dilepton_test("pNb" "3.5")
dilepton_test("pi-p" "0.56")

dilepton_test("CC" "1.0")
dilepton_test("CC" "2.0")

dilepton_test("ArKCl" "1.76")
dilepton_test("AuAu" "1.25")
dilepton_test("AgAg" "1.67")

## tests with filtering ##
if (WITH_EXP_DATA)
  dilepton_test("pp" "1.25" "HadesSingleAcc-p125p-APR06-effgt5-RK-v2.acc" "0.05" "1.8")
  dilepton_test("pp" "2.2"  "HadesSingleAcc-p22p-JAN04-effgt5-RK-v2.acc"  "0.10" "2.0")
  dilepton_test("pp" "3.5"  "HadesSingleAcc-p35p-APR07-effgt5-RK-v2.acc"  "0.08" "2.0")
  dilepton_test("pNb" "3.5"  "HadesSingleAcc-p35Nb-SEP08-effgt5-MVA-v1.acc"  "0.1" "2.0")

  dilepton_test("CC" "1.0" "HadesSingleAcc-c10c-AUG04-effgt5-KP-v2.acc" "0.05" "1.8")
  dilepton_test("CC" "2.0" "HadesSingleAcc-c20c-NOV02-effgt0-KP-v1.acc" "0.05" "1.8")
  dilepton_test("ArKCl" "1.76" "HadesSingleAcc-ar17kcl-SEP05-effge5-HC-flat-v2.acc" "0.1" "1.1")
endif()

## box tests

box_test("omega")
box_test("phi")
box_test("eta-prime")

# BR test

shining_BR_test()


## Top-level target for all dilepton tests that are run regularly within the
## analysis suite. The remaining colliding systems need to be started separately.

add_custom_target(dileptons
                  DEPENDS ${dileptons_ana})
