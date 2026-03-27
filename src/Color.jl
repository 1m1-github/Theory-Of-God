function scalar2rgba(t)
    d = abs(t - ○)
    alpha = T(2) * d          # distance from center = confidence

    # map position to hue angle: full circle over [0,1]
    # t=0 → hue=0 (red), t=0.5 → undefined (doesn't matter, alpha=0)
    # t=1 → hue=180 (cyan), cold/hot sides are opposite on wheel
    hue = mod(t * T(360), T(360))

    # HSV→RGB with S=1, V=1
    c = one(T)
    x = one(T) - abs(mod(hue / T(60), T(2)) - one(T))
    r, g, b = if hue < T(60);      (c, x, zero(T))
    elseif hue < T(120); (x, c, zero(T))
    elseif hue < T(180); (zero(T), c, x)
    elseif hue < T(240); (zero(T), x, c)
    elseif hue < T(300); (x, zero(T), c)
    else;              (c, zero(T), x)
    end

    return (r, g, b, alpha)
end

rgba2scalar(rgba::PNGFiles.ColorTypes.RGBA) = rgba2scalar(rgba.r, rgba.g, rgba.b, rgba.alpha)
function rgba2scalar(r, g, b, a)
    # RGB → hue (standard formula)
    cmax = max(r, g, b)
    cmin = min(r, g, b)
    delta = cmax - cmin

    if iszero(delta)
        # achromatic — use alpha to get distance from 0.5,
        # but we've lost which side. default to 0.5 (no info)
        return ○
    end

    hue = if cmax == r
        T(60) * mod((g - b) / delta, T(6))
    elseif cmax == g
        T(60) * ((b - r) / delta + T(2))
    else
        T(60) * ((r - g) / delta + T(4))
    end

    if hue < zero(T)
        hue += T(360)
    end

    return hue / T(360)
end

a=unique(rgba)
a[3]
a[2]
rgba2scalar(a[3])
rgba2scalar(a[2])
rgba2scalar(0.831,0.831,0.831,1.0)
unique(rgba2scalar.(rgba))