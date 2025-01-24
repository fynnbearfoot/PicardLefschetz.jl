# module ContourIntersection

# using Contour
# using StaticArrays
# using GeometryBasics


export crosses_point, intersection, closest_intersection

struct Line2
    s::StaticArrays.SVector{2,Float64}
    e::StaticArrays.SVector{2,Float64}

    Line2(t1::Tuple{Float64, Float64}, t2::Tuple{Float64, Float64}) =
        new(SA[t1...],SA[t2...])
end



function crosses_point(line::Line2, point::Point, tolerance::Float64=0.25)
    if line.s.x <= line.e.x
        sx = line.s.x
        ex = line.e.x
    else
        ex = line.s.x
        sx = line.e.x
    end

    if line.s.y <= line.e.y
        sy = line.s.y
        ey = line.e.y
    else
        ey = line.s.y
        sy = line.e.y
    end

    sx -= tolerance*abs(sx-ex)
    ex += tolerance*abs(sx-ex)
    sy -= tolerance*abs(sy-ey)
    ey += tolerance*abs(sy-ey)

    xbound = (sx <= point[1] <= ex)
    ybound = (sy <= point[2] <= ey)
    return (xbound && ybound)
end

### maybe these two following functions could be united actually
function crosses_point(c::Curve2, p::Point, tolerance::Float64=0.25)
    crossing = false

    for i in 1:(length(c.vertices)-1)
        l = Line2(c.vertices[i], c.vertices[i+1])
        if crosses_point(l, p , tolerance)
            return i #crossing = true
        end
    end 
    return crossing
end

function closest_intersection(ip::Int64, c::Curve2, p::Point, Δinit::Float64)
    ip_new = deepcopy(ip)
    n = 1
    Δ = Δinit
    while !(ip_new==false) && n < 10
#         @show Δ
        ip = ip_new
        ip_new = crosses_point(c, p, Δ)
#         @show saddle_ip
        n +=1
        Δ -= 0.1
    end
    ip
end



function intersection(l1::Line2, l2::Line2) 
    x1, y1 = l1.s
    x2, y2 = l1.e
    x3, y3 = l2.s
    x4, y4 = l2.e
    
    denom = (y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1)
    
    if denom == 0
        return nothing  # Lines are parallel
    end
    
    ua = ((x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3)) / denom

    ub = ((x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3)) / denom
    
    if 0 <= ua <= 1 && 0 <= ub <= 1
        intersection = (x1 + ua*(x2 - x1), y1 + ua*(y2 - y1))
        return Point(intersection...)  # Lines intersect
    else
        return nothing  # Lines do not intersect
    end

    # a1 = y2 - y1
    # b1 = x1 - x2
    # c1 = a1 * x1 + b1 * y2
 
    # a2 = y4 - y3
    # b2 = x3 - x4
    # c2 = a2 * x3 + b2 * y3
 
    # diff = (y2 - y1) * (x3 - x4) - (y4 - y3) * (x1 - x2)

    # point = Point(
    #     (b2 * c1 - b1 * c2) / diff, 
    #     (a1 * c2 - a2 * c1) / diff)
    # if crosses_point(l1, point) && crosses_point(l2, point)
    #     return point
    # else
    #     return nothing
    # end 
end

function intersection(c1::Curve2, c2::Curve2)
    intersection_points = Vector{Point}()
    for i in 1:(length(c1.vertices)-1), i2 in 1:(length(c2.vertices)-1)
        l1 = Line2(c1.vertices[i],c1.vertices[i+1])
        l2 = Line2(c2.vertices[i2],c2.vertices[i2+1])
        p = intersection(l1,l2)
        if p != nothing
            push!(intersection_points,p)
        end 
         
    end 
    return intersection_points
end 

function intersection(con1_curves::ContourLevel, con2_curves::ContourLevel)
    intersection_points = Vector{Point}()
    for curve1 in lines(con1_curves)
        for curve2 in lines(con2_curves)
            push!(intersection_points,intersection(curve1,curve2)...)
        end 
    end 
    return intersection_points
end 


# end


nothing