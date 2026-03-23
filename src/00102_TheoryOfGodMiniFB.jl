const MAX_RGB = T(mfb_rgb(255, 255, 255))
rgb2c(r, g, b) = T(mfb_rgb(r * 255, g * 255, b * 255)) / MAX_RGB
c2rgb(c2) = begin
    c = floor(UInt32, c2 * MAX_RGB)
    ((c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF) ./ 255
end

window = mfb_open_ex("tog", g.♯..., MiniFB.WF_RESIZABLE)
buffer = zeros(UInt32, prod(g.♯))
MINIFBTASK = @async while true
    yield()
    state = mfb_update(window, buffer)
    state != MiniFB.STATE_OK && break
end
viewer_zero_vs_one_mode = 0
viewer_dim = g.ẑero.d[2]
updatebuffer(g) = begin
    try
    floor.(UInt32, reshape(∃̇(g, Ω), prod(g.♯)) .* MAX_RGB)
    catch e
        showerror(stderr, e, catch_backtrace())
    end
end
const pending_actions = Channel{Function}(32)
# const last_key_time = Ref(0.0)
# const DEBOUNCE_MS = 50
# function keyboard_cb(window::Ptr{Cvoid}, key::Int32, mod::Int32, is_pressed::Bool)::Cvoid
#     try
#         # !is_pressed && return nothing
#         # now = time() * 1000
#         # (now - last_key_time[]) < DEBOUNCE_MS && return nothing
#         # last_key_time[] = now

#         println("Key pressed: $key, $mod, $(is_pressed)")
#         println("viewer_zero_vs_one_mode=$viewer_zero_vs_one_mode")
#         println("viewer_dim=$viewer_dim")
#         println("ẑero.μ=$(g.ẑero.μ)")
#         println("f̂ocus.μ=$(g.f̂ocus.μ)")

#         global viewer_dim
#         global viewer_zero_vs_one_mode

#         if key == 48
#             viewer_zero_vs_one_mode = 1 - viewer_zero_vs_one_mode
#         elseif 49 ≤ key ≤ 49 + 8
#             viewer_dim = g.ẑero.d[key+1-49]
#         elseif key == 265
#             let d = viewer_dim
#                 if iszero(viewer_zero_vs_one_mode)
#                     put!(pending_actions, ĝ -> moveup(ĝ, d))
#                 else
#                     put!(pending_actions, ĝ -> focusup(ĝ, d))
#                 end
#             end
#         elseif key == 264
#             let d = viewer_dim
#                 if iszero(viewer_zero_vs_one_mode)
#                     put!(pending_actions, ĝ -> movedown(ĝ, d))
#                 else
#                     put!(pending_actions, ĝ -> focusdown(ĝ, d))
#                 end
#             end
#         elseif key == 81
#             put!(pending_actions, ĝ -> jerkdown(ĝ))
#         elseif key == 87
#             put!(pending_actions, ĝ -> jerkup(ĝ))
#         elseif key == 69
#             put!(pending_actions, ĝ -> scaledown(ĝ))
#         elseif key == 82
#             put!(pending_actions, ĝ -> scaleup(ĝ))
#         end
#     catch e
#         showerror(stderr, e, catch_backtrace())
#     end
#     global buffer = updatebuffer(g)
#     println("after:")
#     println("viewer_zero_vs_one_mode=$viewer_zero_vs_one_mode")
#     println("viewer_dim=$viewer_dim")
#     println("ẑero.μ=$(g.ẑero.μ)")
#     println("f̂ocus.μ=$(g.f̂ocus.μ)")
# end
function keyboard_cb(window::Ptr{Cvoid}, key::Int32, mod::Int32, is_pressed::Bool)::Cvoid
    try
        # !is_pressed && return
        # now = time() * 1000
        # (now - last_key_time[]) < DEBOUNCE_MS && return
        # last_key_time[] = now

        global g
        global viewer_dim
        global viewer_zero_vs_one_mode
        println("Key pressed: $key, $mod, $(is_pressed)")
        println("viewer_zero_vs_one_mode=$viewer_zero_vs_one_mode")
        println("viewer_dim=$viewer_dim")
        println("ẑero.μ=$(g.ẑero.μ)")
        println("f̂ocus.μ=$(g.f̂ocus.μ)")
        if is_pressed
            if key == 48 # 0
                viewer_zero_vs_one_mode = 1 - viewer_zero_vs_one_mode
            elseif 49 ≤ key ≤ 49 + 8 # numbers
                viewer_dim = g.ẑero.d[key+1-49]
            elseif key == 265 # up
                if iszero(viewer_zero_vs_one_mode)
                    # global g = moveup(g, viewer_dim)
                    # g = moveup(g, viewer_dim)
                    put!(pending_actions, g -> moveup(g, viewer_dim))
                    # put!(pending_actions, ĝ -> moveup(ĝ, viewer_dim))
                else
                    g = focusup(g, viewer_dim)
                end
            elseif key == 264 # down
                if iszero(viewer_zero_vs_one_mode)
                    # global g = movedown(g, viewer_dim)
                    # g = movedown(g, viewer_dim)
                    put!(pending_actions, g -> movedown(g, viewer_dim))
                    # put!(pending_actions, ĝ -> movedown(ĝ, viewer_dim))
                else
                    g = focusdown(g, viewer_dim)
                end
            elseif key == 81 # q
                g = jerkdown(g)
            elseif key == 87 # w
                g = jerkup(g)
            elseif key == 69 # e
                g = scaledown(g)
            elseif key == 82 # r
                g = scaleup(g)
            end
        end
    catch e
        bt = catch_backtrace()
        showerror(stderr, e, bt)
    end
    # global buffer = updatebuffer(g)
    println("after:")
    println("viewer_zero_vs_one_mode=$viewer_zero_vs_one_mode")
    println("viewer_dim=$viewer_dim")
    println("ẑero.μ=$(g.ẑero.μ)")
    println("f̂ocus.μ=$(g.f̂ocus.μ)")
end
kb_cfunc = @cfunction(keyboard_cb, Cvoid, (Ptr{Cvoid}, Int32, Int32, Bool))
ccall((:mfb_set_keyboard_callback, MiniFB.libminifb), Cvoid,
    (Ptr{Cvoid}, Ptr{Cvoid}),
    window, kb_cfunc)
# schedule(MINIFBTASK, InterruptException(), error=true)
# mfb_close(window)
