using Serialization

function params_to_fname()
    # return a string filename for the current param set

    prefix = dimstring *
             "$(nameof(typeof(initial_condition)))" *
             "N$N" *
             "K$K" *
             "VF$(nameof(typeof(volume_flux)))" *
             "BL$blend"

    if blend == :subcell
        prefix *= "KS$(nameof(typeof(knapsack_solver)))"
    end

    middlefix = "TS$(nameof(typeof(timestepper)))" *
                "AD$adaptive" *
                "DT$dt"

    if adaptive
        middlefix *= "AT$abstol" *
                     "RT$reltol"
    end

    postfix = "KSC$(nodewise_shock_capturing)" *
              "YSC$(shock_capturing)"

    return prefix * middlefix * postfix
end

function save_run()
    # Save sol.t and sol.u to file
    if use_run_saver
        serialize("object_transfer/$(params_to_fname())", sol)
    end
end

function load_save_if_exists()
    # Check if the run exists, if so, return sol. Otherwise, return -1
    if use_run_saver && isfile("object_transfer/$(params_to_fname())")
        println("Run with current params exists, it is being loaded into current context.")
        return deserialize("object_transfer/$(params_to_fname())")
    else
        return -1
    end
end