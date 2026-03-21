# @kernel function κ!(rgba, Φ, Φi, ♯)
    # rgba[i] = clamp(0f0, 0f0, 1f0)
    # ix = ((idx - Int32(1)) % W) + Int32(1)
    # iy = ((idx - Int32(1)) ÷ W) + Int32(1)

    # # Uniform X/Y in [0,1]
    # px = Float32(ix - Int32(1)) / Float32(W - Int32(1))
    # py = Float32(iy - Int32(1)) / Float32(H - Int32(1))

    # acc_r = 0f0
    # acc_g = 0f0
    # acc_b = 0f0
    # acc_a = 0f0

    # for iz in Int32(1):D
    #     acc_a >= 0.999f0 && break
    #     pz = grid_z[iz]

    #     for oi in Int32(1):n_objects
    #         b = (oi - Int32(1)) * Int32(FLOATS_PER_OBJ)
    #         shape = scene_data[b+Int32(1)]
    #         cx = scene_data[b+Int32(2)]
    #         cy = scene_data[b+Int32(3)]
    #         cz = scene_data[b+Int32(4)]
    #         p1 = scene_data[b+Int32(5)]
    #         p2 = scene_data[b+Int32(6)]
    #         cr = scene_data[b+Int32(7)]
    #         cg = scene_data[b+Int32(8)]
    #         cb = scene_data[b+Int32(9)]
    #         al = scene_data[b+Int32(10)]
    #         th = scene_data[b+Int32(11)]

    #         lr = 0f0
    #         lg = 0f0
    #         lb = 0f0
    #         la = 0f0

    #         # --- Fusion: 3 separate kernels, one dispatch, compare as Float32 ---
    #         if shape == 1f0
    #             (lr, lg, lb, la) = eval_disc(px, py, pz, cx, cy, cz, p1, th, cr, cg, cb, al)
    #         elseif shape == 2f0
    #             (lr, lg, lb, la) = eval_ring(px, py, pz, cx, cy, cz, p1, p2, th, cr, cg, cb, al)
    #         elseif shape == 3f0
    #             (lr, lg, lb, la) = eval_square(px, py, pz, cx, cy, cz, p1, th, p2, cr, cg, cb, al)
    #         end

    #         if la > 0f0
    #             rem = 1f0 - acc_a
    #             acc_r += lr * la * rem
    #             acc_g += lg * la * rem
    #             acc_b += lb * la * rem
    #             acc_a += la * rem
    #         end
    #     end
    # end

    # # Dark background
    # rem = 1f0 - acc_a
    # image_r[idx] = clamp(acc_r + 0.03f0 * rem, 0f0, 1f0)
    # image_g[idx] = clamp(acc_g + 0.03f0 * rem, 0f0, 1f0)
    # image_b[idx] = clamp(acc_b + 0.06f0 * rem, 0f0, 1f0)
    # image_a[idx] = clamp(acc_a, 0f0, 1f0)
# end

# function render(; W=1920, H=1080, D=256, show_grid=true)
#     println("Rendering $(W)×$(H), D=$(D), all coords in [0,1]³")

#     grid_z = log_z_grid(D)

#     if show_grid
#         println("  Log Z: $(round(grid_z[1],digits=4)) → $(round(grid_z[end],digits=2))")
#         bins = range(0.0, 1.0, length=11)
#         for i in 1:length(bins)-1
#             count = sum(bins[i] .<= grid_z .< bins[i+1])
#             bar = "█"^max(0, round(Int, count * 40 / D))
#             println("    $(rpad("$(round(bins[i],digits=1))-$(round(bins[i+1],digits=1))",7)) $(lpad(count,3)) $bar")
#         end
#     end

#     scene_flat, n_objects = make_scene()
#     gz_gpu = to_gpu(grid_z)
#     scene_gpu = to_gpu(scene_flat)

#     N = W * H
#     # img_r = gpu_zeros(Float32, N)
#     # T=Float32
#     # N=10
#     rgba = KernelAbstractions.zeros(GPU_BACKEND, T, ♯...)
#     # img_g = gpu_zeros(Float32, N)
#     # img_b = gpu_zeros(Float32, N)
#     # img_a = gpu_zeros(Float32, N)
#     render!(GPU_BACKEND, 2^2^3)(
#         rgba,
#         gz_gpu, scene_gpu,
#         Int32(n_objects), Int32(W), Int32(H), Int32(length(grid_z)),
#         ndrange=N,
#     )
#     KernelAbstractions.synchronize(GPU_BACKEND)
#     r = reshape(Array(img_r), W, H)
#     g = reshape(Array(img_g), W, H)
#     b = reshape(Array(img_b), W, H)
#     img = RGB.(permutedims(r, (2, 1)), permutedims(g, (2, 1)), permutedims(b, (2, 1)))
#     return img, t
# end
