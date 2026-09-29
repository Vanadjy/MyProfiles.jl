function f_hists_data(n_probs::Int, n_algos::Int)
    data = fill(Dict(1:n_algos .=> [sort(50*rand(20), rev = true) for i in 1:n_algos]), n_probs)
    return data
end

function N_hists_data(n_probs::Int, n_algos::Int)
    data = fill(Dict(1:n_algos .=> [sort(rand(1:500, 20)) for i in 1:n_algos]), n_probs)
    return data
end