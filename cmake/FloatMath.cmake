function(float_math VAR expr)
    execute_process(COMMAND "python3" "-Sc" "print(${expr})"
                    OUTPUT_VARIABLE output
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(${VAR} ${output} PARENT_SCOPE)
endfunction()

