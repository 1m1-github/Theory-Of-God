# const MINIFB_WINDOW = Ref(mfb_open_ex("tog", G[].♯..., MiniFB.WF_RESIZABLE))
# const MINIFB_BUFFER = Ref(zeros(UInt32, prod(G[].♯)))
# const MINIFB_TASK = Threads.@spawn while true
#     yield()
#     state = mfb_update(MINIFB_WINDOW[], MINIFB_BUFFER[])
#     state != MiniFB.STATE_OK && break
# end
# const PENDING_ACTIONS = Channel{Function}(32)
# const CHANGE_MODE = Ref(0) # 0=zero, 1=focus, 2=ρ
# const CHANGE_DIM_INDEX = Ref(2)
const MAX_RGB = T(mfb_rgb(255, 255, 255))
rgb2c(r, g, b) = T(mfb_rgb(r * 255, g * 255, b * 255)) / MAX_RGB
rgb2c(r, g, b, a) = rgb2c(r, g, b)
rgb2c(rgba) = rgb2c(rgba.r, rgba.g, rgba.b)
c2rgb(c) = begin
    c = floor(UInt32, c * MAX_RGB)
    ((c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF) ./ 255
end
# function keyboard_cb(window::Ptr{Cvoid}, key::Int32, mod::Int32, is_pressed::Bool)::Cvoid
#     try # todo rm
#         println("Key pressed: $key, $mod, $(is_pressed)")
#         println("CHANGE_MODE=$CHANGE_MODE[]")
#         println("CHANGE_DIM_INDEX=$CHANGE_DIM_INDEX[]")
#         println("ẑero.μ=$(G[].ẑero.μ)")
#         println("ône.μ=$(G[].ône.μ)")
#         println("ône.μ=$(G[].ρ)")
#         if key == 48
#             global CHANGE_MODE[] = (CHANGE_MODE[] + 1) % 3
#         elseif 49 ≤ key ≤ 49 + 8
#             global CHANGE_DIM_INDEX[] = key - 48
#         elseif key == 265
#             if CHANGE_MODE[] == 0
#                 put!(PENDING_ACTIONS, g -> moveup(g, CHANGE_DIM_INDEX[]))
#             elseif CHANGE_MODE[] == 1
#                 put!(PENDING_ACTIONS, g -> focusup(g, CHANGE_DIM_INDEX[]))
#             else
#                 put!(PENDING_ACTIONS, g -> scaleup(g, CHANGE_DIM_INDEX[]))
#             end
#         elseif key == 264
#             if CHANGE_MODE[] == 0
#                 put!(PENDING_ACTIONS, g -> movedown(g, CHANGE_DIM_INDEX[]))
#             elseif CHANGE_MODE[] == 1
#                 put!(PENDING_ACTIONS, g -> focusdown(g, CHANGE_DIM_INDEX[]))
#             else
#                 put!(PENDING_ACTIONS, g -> scaledown(g, CHANGE_DIM_INDEX[]))
#             end
#         elseif key == 81
#             put!(PENDING_ACTIONS, g -> jerkdown(g))
#         elseif key == 87
#             put!(PENDING_ACTIONS, g -> jerkup(g))
#         end
#     catch e
#         showerror(stderr, e, catch_backtrace())
#     end
#     println("after:")
#     println("CHANGE_MODE=$CHANGE_MODE[]")
#     println("CHANGE_DIM_INDEX=$CHANGE_DIM_INDEX[]")
#     println("ẑero.μ=$(G[].ẑero.μ)")
#     println("ône.μ=$(G[].ône.μ)")
#     println("ône.μ=$(G[].ρ)")
# end
# kb_cfunc = @cfunction(keyboard_cb, Cvoid, (Ptr{Cvoid}, Int32, Int32, Bool))
# ccall((:mfb_set_keyboard_callback, MiniFB.libminifb), Cvoid,
#     (Ptr{Cvoid}, Ptr{Cvoid}),
#     MINIFB_WINDOW[], kb_cfunc)
# schedule(MINIFBTASK, InterruptException(), error=true)
# mfb_close(MINIFB_WINDOW[])


# updatebuffer() = begin
#     try # todo rm
#         Base.invokelatest() do
#             global MINIFB_BUFFER[] = floor.(UInt32, reshape(∃̇(G[], Ω[]), prod(G[].♯)) .* MAX_RGB)
#         end
#     catch e
#         showerror(stderr, e, catch_backtrace())
#     end
# end
# const UPDATE_MINIFB_BUFFER_TASK = Threads.@spawn while true
#     yield()
#     # sleep(1) # todo rm
#     updatebuffer()
# end

# const CHANGE_TASK = Threads.@spawn while true
#     yield()
#     while isready(PENDING_ACTIONS)
#         Α = take!(PENDING_ACTIONS)
#         global G[] = Α(G[])
#     end
# end