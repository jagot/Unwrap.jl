using PyCall
pygui(:qt)
using PyPlot
matplotlib[:rcdefaults]()
ion()

using Unwrap
using DSP

x = linspace(0,10,3000)*2π
φ = angle(exp(im*(x-8π).^2))
φφ = @time unwrap1d(φ)
φφφ = @time unwrap(φ)

figure("1d",figsize=(5,4))
clf()

subplot(211)
plot(x/2π,φ/π)
margins(0,0.1)
gca()[:set_xticklabels]([])
ylabel(L"Wrapped [$\pi$]")
subplot(212)
plot(x/2π,φφ/π,label="Unwrap.jl")
plot(x/2π,φφφ/π,label="DSP.jl")
margins(0,0.1)
ylabel(L"Unwrapped [$\pi$]")
legend(framealpha=0.75)
tight_layout()
savefig("1d.png",dpi=72,transparent=true)
