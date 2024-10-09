using DataFrames, Random

function get_evts_index(df::DataFrame)
    event_number = df[1,:evt]
    evt_counter = 1
    evts_info = []
    start_evt_index = 1
    for i = 1 : 1 : length(df[!,:evt])
        if  event_number == df[i,:evt]
            continue
        else
            end_evt_index = i - 1
            event_number = df[i,:evt]
            push!(evts_info,[evt_counter start_evt_index end_evt_index])
            evt_counter += 1
            start_evt_index = i
        end
    end
    #evt_counter += 1
    start_evt_index = end_evt_index + 1
    end_evt_index = length(df[!,:evt])
    push!(evts_info,[evt_counter start_evt_index end_evt_index])
    return vcat(evts_info...)
end

function greet_DUNEatLapp()
    return "Hello DUNEatLapp!"
end


function test_function()
	return "Haut les mains peau de lapin !"
end

