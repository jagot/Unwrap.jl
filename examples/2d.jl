using PyCall
pygui(:qt)
using PyPlot
matplotlib[:rcdefaults]()
ion()

using Unwrap

x = linspace(-1,1,200)
y = x
r² = broadcast(+, x.^2, (y').^2)
z = exp(-r²)
φ = angle(exp(im*z*10π))

φ_old = copy(φ)
φ1 = unwrap2d(φ)
φ2 = unwrap2d(φ1) # Iterating improves result

function plot_map(args...;kwargs...)
    p = pcolormesh(args...;
                   rasterized=true,
                   cmap=plt[:cm][:get_cmap]("viridis"),
                   kwargs...)
    margins(0,0)
    p
end

function plot_phase(φ,gs,i,label,ylabels=false)
    fig[:add_subplot](get(gs, (1,i-1)))
    pm = plot_map(x,y,φ/2π)
    xlabel(L"x")
    if ylabels
        ylabel(L"y")
    else
        gca()[:set_yticklabels]([])
    end
    cax=fig[:add_subplot](get(gs, (0,i-1)))
    colorbar(pm,cax=cax,orientation="horizontal",label=label)
    cax[:set_xticklabels](cax[:get_xticklabels](),rotation=45)
    cax[:xaxis][:tick_top]()
    cax[:xaxis][:set_ticks_position]("both")
    cax[:xaxis][:set_label_position]("top")
end
fig = figure("2d",figsize=(10,4.5))
clf()
gs = matplotlib[:gridspec][:GridSpec](2,3,height_ratios=(1,10))
plot_phase(φ_old,gs,1,L"Wrapped [$\pi$]",true)
plot_phase(φ1,gs,2,L"1st unwrapped [$\pi$]")
plot_phase(φ2,gs,3,L"2nd unwrapped [$\pi$]")
gs[:tight_layout](fig)

savefig("2d.png",dpi=72,transparent=true)
