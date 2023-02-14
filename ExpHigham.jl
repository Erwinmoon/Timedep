using LinearAlgebra

## run code by include("ExpHigham.jl")

## A = Array{BigFloat} ; A = rand(BigFloat , (3 , 3));

function exph(A::Array{Float64})::Array{Float64}
    ## culculate exp(A) by higham methods, 
    ## where A is a n by n matrix
    (s , m , A2 , A4 , A6) = exphparams(A);
    A = A / (2^s);
    A2 = A2 / (2 ^ (2*s));
    A4 = A4 / (2 ^ (4*s));
    A6 = A2 / (2 ^ (6*s));
    
    F = padeApprox(A, A2, A4, A6, m);

    for k in 1:s
        F = F * F;
    end
    return F
end

function exphparams(A::Array{Float64})
    ## expm_params Obtain scaling parameter and order of the Pade approximant.
    T = Float64;
    ## Coefficients of backwards error function.
    coeff = Vector{T}(undef, 5);
    coeff[1] = 1/100800;
    coeff[2] = 1/10059033600;
    coeff[3] = 1/4487938430976000;
    coeff[4] = 1/5914384781877411840000;
    coeff[5] = 1/113250775606021113483283660800000000; 

    s = 0; ## rescale number

    theta = Vector{T}(undef, 5);
    ##  m_val is one of [3 5 7 9 13];
    ## theta_m for m=1:13.

    theta = [## 3.650024139523051e-008
    ## 5.317232856892575e-004
    1.4955852179582920e-002;  ## m_vals = 3
    ## 8.536352760102745e-002
    2.5393983300632300e-001;  ## m_vals = 5
    ## 5.414660951208968e-001
    9.5041789961629320e-001;  ## m_vals = 7
    ## 1.473163964234804e+000
    2.0978479612570680e+000;  ## m_vals = 9
    ## 2.811644121620263e+000
    ## 3.602330066265032e+000
    ## 4.458935413036850e+000
    5.3719203511481520e+000];  ## m_vals = 13

    A2 = A * A;
    A4 = A2 * A2;
    A6 = A2 * A4;
    d4 = opnorm(A4 , 1)^(1 / 4);
    d6 = opnorm(A6 , 1)^(1 / 6);
    d8 = opnorm(A4 * A4 , 1)^(1 / 8);
    d10 = opnorm(A4 * A6 , 1)^(1 / 10);
    eta1 = max(d4 , d6);
    eta3 = max(d6 , d8);
    eta4 = max(d8 , d10);
    eta5 = min(eta3 , eta4);
    s = max(cld(log2(eta5 / theta[5]) , 1) , 0);
    s = s + ell(A ./ (2^s) , coeff[5] , 13);
    m = 13;

    return s, m, A2, A4, A6
end

function ell(A::Array{Float64}, coeff::Float64, m_val::Int)::Int
    scaledA = coeff .^ (1 / (2 * m_val + 1)) .* abs.(A);
    alpha = norm(scaledA , 2 * m_val + 1) / opnorm(A , 1);
    t = max(cld(log2(2 * alpha / eps(typeof(alpha))) / (2 * m_val) , 1) , 0);
    return t
end

function getPadeCoefficients(m::Int)::Array{Float64}
    ## get_pade_coefficients Coefficients of numerator P of Pade approximant
    ## C = get_pade_coefficients returns coefficients of numerator
    ## of [m/m] Pade approximant, where m = 3,5,7,9,13.
    if m == 3
        c = [120, 60, 12, 1];
    elseif m == 5
        c = [30240, 15120, 3360, 420, 30, 1];
    elseif m == 7
        c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
    elseif m == 9
        c = [17643225600, 8821612800, 2075673600, 302702400, 30270240,
        2162160, 110880, 3960, 90, 1];
    else  ## m == 13
        c = [64764752532480000, 32382376266240000, 7771770303897600, 
            1187353796428800,  129060195264000,   10559470521600, 
            670442572800,      33522128640,       1323241920,
            40840800,          960960,            16380,  182,  1];
    end
    return c
end

function padeApprox(A::Array{Float64} , A2::Array{Float64} , A4::Array{Float64} , A6::Array{Float64} , m::Int)::Array{Float64}
    ## pade_approx Computes the Pade approximant to exp(T) of order [m/m].  
    c = getPadeCoefficients(m);
    (n , ~) = size(A);
    E = Matrix{Float64}(I, n, n);
    U = A * (A6*(c[14]*A6 + c[12]*A4 + c[10]*A2) + c[8]*A6 + c[6]*A4 + 
        c[4]*A2 + c[2]*E);
    V = A6*(c[13]*A6 + c[11]*A4 + c[9]*A2) + c[7]*A6 + c[5]*A4 + 
        c[3]*A2 + c[1]*E;  
    F = (V - U) \ (2 * U) + E;
    return F      
end

## main function 
A = Array{Float64};
A = 0.001 * rand(100 , 100);
errA = (exph(A) - exp(A)) / opnorm(exp(A) , 1);
err = opnorm(errA , 1);
println("relative error:  ", err)