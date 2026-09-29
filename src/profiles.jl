export conv_plot, trajectory_plot, f_star, f_0, accuracy, Nap, rap, perf_profile!, data_profile!, accuracy_profile!

function conv_plot(f_hists, N_hists, prob::Int; logscale::Bool = false)
    graph = plot()
    ns = length(keys(N_hists[prob]))
    @inbounds for algo in keys(N_hists[prob])
        if logscale
            plot!(N_hists[prob][algo], f[prob][algo], linetype=:steppre, xaxis=:log10, yaxis=:log10, label=key(N_hists[prob])[i],
                  xlabel="Number of F evaluations",
                  ylabel="F value")
        else
            plot!(N[i], f[i](p), linetype=:steppre, label="algorithm $i")
            xlabel!("Number of F evaluations")
            ylabel!("F value")
        end
        title!("convergence plot for problem $p")
    end
    return graph
end

function trajectory_plot(f, f_hist, xk_hist, yk_hist; 
                        UB = 10.0, 
                        LB = -10.0, 
                        xlabel = true, 
                        ylabel = true,
                        color_bar = true,
                        color_bar_title = true,
                        nb_of_levels::Int = 30,
                        granul::Int = 100,
                        label::String = "Algo",
                        linewidth::Real = 2,
                        color = :red)
    # Trajectory plots available only for dimension 2 or 1.
    @assert length(xk_hist[1]) == 1 "Error: Trajectory plots available only for dimension 2 or 1"
    @assert length(yk_hist[1]) <= 1 "Error: Trajectory plots available only for dimension 2 or 1"
    @assert length(f_hist) == length(xk_hist) "Error: Dimension mismatch ; historics must have same lengths"
    @assert length(f_hist) == length(yk_hist) "Error: Dimension mismatch ; historics must have same lengths"

    # Plots contours of function f
    X = collect(range(LB, UB, granul))
    Y = collect(range(LB, UB, granul))
    plt = contour(X, Y, f;
              levels = nb_of_levels,
              fill = true,
              color = :viridis,
              colorbar = color_bar,
              xlabel = xlabel ? "x" : "",
              ylabel = ylabel ? "y" : "",
              colorbar_title = color_bar_title ? "Levels of function f" : ""
    )
    # Trajectory
    plot!(plt, xk_hist', yk_hist'; color = color, linewidth = linewidth, label = label)
    scatter!(plt, xk_hist', yk_hist'; color = color, markersize=4, markershape=:circle, label="")
    # Starting and end points
    scatter!(plt, [xk_hist[1]], [yk_hist[1]];
            label="Start", color=:red, markersize=6, markershape=:star5)
    scatter!(plt, [xk_hist[end]], [yk_hist[end]];
            label="End", color=:green, markersize=6, markershape=:star5)
    return plt
end

function f_star(f_hists, prob::Int, algo_list::Union{Vector{Int}, Vector{String}})
    n_algos = length(algo_list)

    @assert n_algos > 1 "Trying to compare only one algorithm for performace/data profiles"

    best_val = length(f_hists[prob][algo_list[1]]) == 0 ? Inf : f_hists[prob][algo_list[1]][length(f_hists[prob][algo_list[1]])]
    for a in 2:n_algos
        candidate_val = length(f_hists[prob][algo_list[a]]) == 0 ? Inf : f_hists[prob][algo_list[a]][length(f_hists[prob][algo_list[a]])]
        if best_val > candidate_val
            best_val = candidate_val
        end
    end
    return best_val
end

function f_0(f_hists, prob::Int, algo::Union{Int, String})
    f0 = length(f_hists[prob][algo]) == 0 ? Inf : f_hists[prob][algo][1]
    return f0
end

function f_0(f_hists, prob::Int, algo_list::Union{Vector{Int}, Vector{String}})
    n_algos = length(algo_list)
    @assert n_algos > 1 "Trying to compare only one algorithm for performace/data profiles"

    f0 = length(f_hists[prob][algo_list[1]]) == 0 ? -Inf : f_hists[prob][algo_list[1]][1]
    for a in 2:n_algos
        # Routine to select the highest first feasible f0
        candidate_f0 = length(f_hists[prob][algo_list[a]]) == 0 ? -Inf : f_hists[prob][algo_list[a]][1]
        if f0 < candidate_f0
            f0 = candidate_f0
        end
    end
    # If f0 is -Inf, then no algorithm found a feasible point, so we set f0 to NaN
    f0 = (f0 == -Inf ? NaN : f0)
    return f0
end


function accuracy(f_hists, k::Int, prob::Union{Int, String}, algo::Union{Int, String}, algo_list::Union{Vector{Int}, Vector{String}})
    f_N = length(f_hists[prob][algo]) == 0 ? Inf : f_hists[prob][algo][k]
    return ((f_N - f_0(f_hists, prob, algo_list))/(f_star(f_hists, prob, algo_list) - f_0(f_hists, prob, algo_list)))
end

function accuracy(f_hists, k::Int, prob::Union{Int, String}, algo::Union{Int, String}, algo_list::Union{Vector{Int}, Vector{String}}, f_star::Real)
    f_N = length(f_hists[prob][algo]) == 0 ? Inf : f_hists[prob][algo][k]
    return ((f_N - f_0(f_hists, prob, algo_list))/(f_star - f_0(f_hists, prob, algo_list)))
end

function Nap(f_hists, N_hists, algo::Union{Int, String}, prob::Union{Int, String}, τ::Real, algo_list::Union{Vector{Int}, Vector{String}})
    Nap = Inf
    Tap = false
    i = 1
    @assert length(N_hists[prob][algo]) == length(f_hists[prob][algo]) "ERROR: The length of N_hist and f_hist should be the same for problem $prob and algorithm $algo."
    while !(Tap || i >= length(N_hists[prob][algo]))
        i += 1
        if accuracy(f_hists, i, prob, algo, algo_list) ≥ 1 - τ
            Nap = N_hists[prob][algo][i]
            Tap = true
        end
    end
    return Nap, Tap
end

function rap(f_hists, N_hists, algo::Union{Int, String}, prob::Union{Int, String}, τ::Real, algo_list::Union{Vector{Int}, Vector{String}})
    Nap_ref, Tap_ref = Nap(f_hists, N_hists, algo, prob, τ, algo_list)
    champ_Nap = Nap_ref
    rap = Inf
    if Tap_ref
        for alg in algo_list
            Nap_algo, Tap_algo = Nap(f_hists, N_hists, alg, prob, τ, algo_list)
            if Tap_algo && (Nap_algo < champ_Nap)
                champ_Nap = Nap_algo
            end
        end
        rap = Nap_ref / champ_Nap
    end
    return rap
end

function perf_profile!(y, αs, f_hist, N_hist, prob_list::Vector{Int}, algo::Union{Int, String}, τ::Real, algo_list::Union{Vector{Int}, Vector{String}})
    count = 0
    @inbounds for l in eachindex(αs)
        α = αs[l]
        @inbounds for prob in eachindex(prob_list)
            if (rap(f_hist, N_hist, algo, prob, τ, algo_list) ≤ α)
                count += 1
            end
        end
        ρ = count / (length(prob_list))
        y[l] = ρ
        count = 0
    end
    return y
end

function data_profile!(y, ks, f_hist, N_hist, prob_list::Vector{Int}, algo::Union{Int, String}, τ::Real, algo_list::Union{Vector{Int}, Vector{String}}; λ_toggle::Bool = false, effort_choice::String = "UL")
    count = 0
    @inbounds for l in eachindex(ks)
        k = ks[l]
        @inbounds for prob in eachindex(prob_list)
            Nap_data, Tap_data = Nap(f_hist, N_hist, algo, prob, τ, algo_list)
            model = get_bilevel_problem(prob_list[prob])
            if λ_toggle # if we scaled the UL evaluations with λ
                if effort_choice == "UL"
                    dimprob = (model.dim[1] + 1)
                elseif effort_choice == "LL"
                    dimprob = (model.dim[2] + 1)
                else
                    dimprob = (model.dim[1] + 1) * (model.dim[2] + 1)
                end
            else # otherwise, depends on the budget choice
                if effort_choice == "UL"
                    dimprob = model.dim[1] + 1
                elseif effort_choice == "LL"
                    dimprob = model.dim[2] +1
                else
                    dimprob = (model.dim[1] + 1) * (model.dim[2] + 1)
                end
            end
            if Nap_data ≤ k * (dimprob) * Tap_data
                count += 1
            end
        end
        dk = count / (length(prob_list))
        y[l] = dk
        count = 0
    end
    return y
end

function accuracy_profile!(y, ds, f_hist, prob_list::Vector{Int}, algo::Union{Int, String}, algo_list::Union{Vector{Int}, Vector{String}}; opt_known::Bool = false)
    count = 0
    @inbounds for i in eachindex(ds)
        d = ds[i]
        @inbounds for prob in eachindex(prob_list)
            k = length(f_hist[prob][algo])
            if opt_known
                bilevel_prob = get_bilevel_problem(prob_list[prob])
                f_star = get_opt_val(bilevel_prob)
                f_acc_tot = accuracy(f_hist, k, prob, algo, algo_list, f_star)
            else
                f_acc_tot = accuracy(f_hist, k, prob, algo, algo_list)
            end
            if isinf(f_acc_tot) || isnan(f_acc_tot)
                f_acc_tot = 0.0
            end
            if  -log10(1 - f_acc_tot) ≥ d
                count += 1
            end
        end
        ratio = count / length(prob_list)
        y[i] = ratio
        count = 0
    end
    return y
end