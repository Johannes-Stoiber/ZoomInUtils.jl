module ZoomInUtils

using Healpix
using GadgetIO
using LinearAlgebra

function build_ZoomIn_boundary_map(snap_base::String, center::Vector; NSIDE::Integer=64)
    # buffer in kpc/h (i.e. code units)
    println("Building Boundary maximal radius map for NSIDE $NSIDE")
    res = Resolution(NSIDE) 
    n_pixels = nside2npix(NSIDE)
    boundary_map = fill(Inf, n_pixels) 
    
    pos_boundary = read_block(snap_base, "POS", parttype = 2)
    for i in 1:size(pos_boundary, 2)
        dx = pos_boundary[1,i] - center[1] 
        dy = pos_boundary[2,i] - center[2] 
        dz = pos_boundary[3,i] - center[3]
        r = sqrt(dx^2 + dy^2 + dz^2)

        pix = vec2pixRing(res, Float64(dx), Float64(dy), Float64(dz))
        if r < boundary_map[pix] 
            boundary_map[pix] = r  
        end     
    end
    return boundary_map 
end 

function filter_ZoomIn_boundary_map!(p::Dict, snap::String, center::Vector, boundary_map::Vector; buffer::Float64=1000.)
    NSIDE = Int64(sqrt(size(boundary_map, 1)/12))
    res = Resolution(NSIDE) 

    pos = read_block(snap, "POS", parttype = 0)
    r = norm.(eachcol( pos .- center )) 
    N = size(r,1)
    valid = trues(N)
    hp = [ vec2pixRing(res, pos[1,i], pos[2,i], pos[3,i]) for i in 1:size(pos,2) ] 
    valid = [r[i] < boundary_map[hp[i]] - buffer for i in 1:size(pos,2) ]
    for key in keys(p)
        p[key] = p[key][valid]
    end 
end 

export build_ZoomIn_boundary_map
export filter_ZoomIn_boundary_map!

end # module ZoomInUtils
