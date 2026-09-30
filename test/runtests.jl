using MyProfiles
using Test

include("create_data.jl")

n_probs = 50
n_algos = 2
prob_list = collect(1:n_probs)
dim_prob_list = fill((2,2), n_probs)
algo_list = collect(1:n_algos)
data_f = f_hists_data(n_probs, n_algos)
data_N = N_hists_data(n_probs, n_algos)

τ = 1e-2

@testset "Computing accuracy and Nap" begin
  for prob in prob_list
    f_star_all = f_star(data_f, prob, algo_list)
    f_0_all = f_0(data_f, prob, algo_list)

    # Tests relative to computing accuracy
    @test f_star_all <= f_0_all
    for algo in algo_list
      f_0_unique = f_0(data_f, prob, algo)
      for k in 1:length(data_f[prob][algo])
        accuracy_test = accuracy(data_f, k, prob, algo, algo_list)

        @test f_0_unique <= f_0_all
        @test f_star_all <= f_0_unique
        @test accuracy_test >= 0
        @test accuracy_test <= 1
      end

      #Tests relative to Nap and rap
      Nap_data, Tap_data = Nap(data_f, data_N, algo, prob, τ, algo_list)
      rap_data = rap(data_f, data_N, algo, prob, τ, algo_list)

      @test Tap_data ? rap_data >= 1 : rap_data == Inf
    end
  end
end

@testset "Computing profiles" begin
  ks = collect(0:1:500)
  y_data = zeros(Float64, length(ks))

  for algo in algo_list
    data_profile!(y_data, 
                  ks, 
                  data_f, 
                  data_N, 
                  prob_list, 
                  dim_prob_list, 
                  algo, 
                  τ, 
                  algo_list; 
                  effort_choice = "Agregate")
    
  end
  for y in y_data
    @test y <= 1.0
  end
end
