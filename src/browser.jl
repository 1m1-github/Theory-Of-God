struct BatchProcessor{T}
    pending::Channel{T}
    notify::Channel{Nothing}
    BatchProcessor{T}() where T = new(Channel{T}(Inf), Channel{Nothing}(1))
end
import Base.put!
function put!(bp::BatchProcessor{T}, item::T) where T
    put!(bp.pending, item)
    isready(bp.notify) || put!(bp.notify, nothing)
end
function start!(f, bp::BatchProcessor{T}) where T
    while true
        take!(bp.notify)
        while true
            batch = T[]
            while isready(bp.pending)
                push!(batch, take!(bp.pending))
            end
            isempty(batch) && break
            # todo add attention
            f(batch)
        end
    end
end

"Serve and execute JavaScript on an HTTP client using SSE"
struct BroadcastBrowser
    stream::HTTP.Streams.Stream
    width::Int
    height::Int
    processor::BatchProcessor{String}
    BroadcastBrowser(stream, width, height) = new(stream, width, height, BatchProcessor{String}())
end
BROADCASTBROWSER = nothing
"`put!(BroadcastBrowser, js)` runs the js on all connected browsers"
put!(::Type{BroadcastBrowser}, js) = put!(BROADCASTBROWSER.processor, js)

const HTML = raw"""
<!DOCTYPE html>
<html>
<body>
<script>
const sse = new EventSource(`/events?width=${document.documentElement.clientWidth}&height=${document.documentElement.clientHeight}`)
sse.onmessage = (e) => eval(e.data)
document.addEventListener('keydown', (e) => {
    fetch('/keypress', {
        method: 'POST',
        body: e.key
    })
})
</script>
</body>
</html>
"""

function safe_write(stream, js)
    try
        HTTP.write(stream, js)
        flush(stream)
        true
    catch e
        e isa Base.IOError || rethrow()
        false
    end
end

function handle_sse(a)
    HTTP.setstatus(a.stream, 200)
    HTTP.setheader(a.stream, "Content-Type" => "text/event-stream")
    HTTP.setheader(a.stream, "Cache-Control" => "no-cache")
    HTTP.startwrite(a.stream)
    start!(a.processor) do input
        for js = input
            js = replace(js, "\n" => ";")
            safe_write(a.stream, "data: $js\n\n") || return
        end
    end
end

function freeport(hint)
    port, server = listenany(hint)
    close(server)
    Int(port)
end

function start(root::Function, keypress::Function, port=freeport(8888))
    HTTP.serve("0.0.0.0", port; stream=true) do stream
        uri = URI(stream.message.target)
        if uri.path == "/"
            HTTP.setstatus(stream, 200)
            HTTP.setheader(stream, "Content-Type" => "text/html")
            HTTP.startwrite(stream)
            HTTP.write(stream, HTML)
        elseif uri.path == "/events"
            params = queryparams(uri)
            width = parse(Int, params["width"])
            height = parse(Int, params["height"])
            @show width, height
            global BROADCASTBROWSER = BroadcastBrowser(stream, width, height)
            root(BROADCASTBROWSER)
            handle_sse(BROADCASTBROWSER)
        elseif uri.path == "/keypress"
            keypress(String(read(stream)))
            HTTP.setstatus(stream, 204)
            HTTP.startwrite(stream)
        else
            HTTP.setstatus(stream, 404)
            HTTP.startwrite(stream)
        end
    end
end
