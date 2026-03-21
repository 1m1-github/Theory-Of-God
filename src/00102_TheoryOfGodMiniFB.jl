using MiniFB

window = mfb_open_ex("tog", g.♯..., MiniFB.WF_RESIZABLE)
buffer = zeros(UInt32, prod(g.♯))
MINIFBTASK = @async while true
    yield()
    state = mfb_update(window, buffer)
    state != MiniFB.STATE_OK && break
end
viewer_zero_vs_one_mode = 0
viewer_dim = g.ẑero.d[2]
function keyboard_cb(window::Ptr{Cvoid}, key::Int32, mod::Int32, is_pressed::Bool)::Cvoid
    try
        global g
        global viewer_dim
        global viewer_zero_vs_one_mode
        # global buffer
        println("Key pressed: $key, $mod, $(is_pressed)")
        println("viewer_zero_vs_one_mode=$viewer_zero_vs_one_mode")
        println("viewer_dim=$viewer_dim")
        println("ẑero.μ=$(g.ẑero.μ)")
        println("ône.μ=$(g.ône.μ)")
        if is_pressed
            if key == 48
                viewer_zero_vs_one_mode = 1 - viewer_zero_vs_one_mode
            elseif 49 ≤ key ≤ 49 + 8
                viewer_dim = g.ẑero.d[key+1-49]
            elseif key == 265
                if iszero(viewer_zero_vs_one_mode)
                    g = moveup(g, viewer_dim)
                else
                    g = focusup(g, viewer_dim)
                end
                # p♯ = calc_p♯(g)
                # buffer = floor.(UInt32, reshape(p♯, prod(g.♯)) .* MAX_RGB)
                # global buffer = c2rgb.(buffer)
                # global buffer = ifelse.(buffer .> 0, 0xFFFFFFFF, 0x00000000)
            elseif key == 264
                if iszero(viewer_zero_vs_one_mode)
                    g = movedown(g, viewer_dim)
                else
                    g = focusdown(g, viewer_dim)
                end
                # p♯ = calc_p♯(g)
                # buffer = floor.(UInt32, reshape(p♯, prod(g.♯)) .* MAX_RGB)
                # global buffer = c2rgb.(buffer)
                # global buffer = ifelse.(buffer .> 0, 0xFFFFFFFF, 0x00000000)
            elseif key == 81
                g = jerkdown(g)
            elseif key == 87
                g = jerkup(g)
            elseif key == 69
                g = scaledown(g)
            elseif key == 82
                g = scaleup(g)
            end
        end
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
        println("after:")
        println("viewer_zero_vs_one_mode=$viewer_zero_vs_one_mode")
        println("viewer_dim=$viewer_dim")
        println("ẑero.μ=$(g.ẑero.μ)")
        println("ône.μ=$(g.ône.μ)")
    end
    return nothing
end
kb_cfunc = @cfunction(keyboard_cb, Cvoid, (Ptr{Cvoid}, Int32, Int32, Bool))
ccall((:mfb_set_keyboard_callback, MiniFB.libminifb), Cvoid,
    (Ptr{Cvoid}, Ptr{Cvoid}),
    window, kb_cfunc)
# schedule(MINIFBTASK, InterruptException(), error=true)
# mfb_close(window)


# using MiniFB

# const WIDTH = 800
# const HEIGHT = 600

# buffer = zeros(UInt32, WIDTH * HEIGHT)

# window = mfb_open("Key Press Viewer", WIDTH, HEIGHT)

# function keyboard_cb(window::Ptr{Cvoid}, key::Int32, mod::Int32, is_pressed::Int32)::Cvoid
#     if is_pressed != 0
#         println("Key pressed: key=$key  mod=$mod")
#     else
#         println("Key released: key=$key  mod=$mod")
#     end
#     return nothing
# end

# const kb_cfun = @cfunction(keyboard_cb, Cvoid, (Ptr{Cvoid}, Int32, Int32, Int32))

# mfb_set_keyboard_callback(window, kb_cfun)

# while mfb_wait_sync(window)
#     state = mfb_update(window, buffer)
#     if state != MiniFB.STATE_OK
#         break
#     end
# end