using LinearAlgebra

## run code by include("ExpTrefethen.jl")

## A = Array{BigFloat} ; A = rand(BigFloat , (3 , 3));

function expt(A::Array{Float64} , N::Int)::Array{Float64} 
    theta = pi * (1 : 2 : N - 1) / N; ## quad pts in (0,pi)

    ## contour can ba change, i.e. Parabolic contour be used here
    z = N * (0.1309 .- 0.1194 * theta.^2 + 0.2500im * theta);
    w = N * (-0.1194 * 2 * theta .+ 0.2500im);
    c = (1im / N) * exp.(z) .* w;

    expA = zero(A); ## value of exp(A)
    (n , m) = size(A);
    E = Matrix{Float64}(I, n, n);

    for k in 1 : 1 : Int(N/2)
        expA = expA - c[k] .* inv((z[k] * E - A));
    end

    expA = 2 * real(expA);
    return expA
end
function exptbig(A::Array{BigFloat} , N::Int)::Array{Float64} 
    theta = pi * (1 : 2 : N - 1) / N; ## quad pts in (0,pi)

    ## contour can ba change, i.e. Parabolic contour be used here
    z = N * (0.1309 .- 0.1194 * theta.^2 + 0.2500im * theta);
    w = N * (-0.1194 * 2 * theta .+ 0.2500im);
    c = (1im / N) * exp.(z) .* w;

    expA = zero(A); ## value of exp(A)
    (n , m) = size(A);
    E = Matrix{BigFloat}(I, n, n);

    for k in 1 : 1 : Int(N/2)
        expA = expA - c[k] .* inv((z[k] * E - A));
    end

    expA = 2 * real(expA);
    return expA
end

## main function 
A = Array{BigFloat};
A = 100 * rand(10 , 10);
B = convert(Array{BigFloat}, A);
N = 64;
errA = (exptbig(B , N) - expt(A , N)) / opnorm(exptbig(B , N) , 1); ## with trefethen
## errA = (exptbig(B , N) - exph(A)) / opnorm(exptbig(B , N) , 1); ## with higham
err = opnorm(errA , Inf);
println("relative error:  ", err)