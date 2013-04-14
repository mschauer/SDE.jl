
function subleading(a, n)
 d = map(D -> 1:D, size(a)[1:end-1])
 sub(a, d..., n)

end

function subleading(a::Matrix, n)
 sub(a, 1:size(a,1), n)
end


function leading(a, n)
 d = map(D -> 1:D, size(a)[1:end-1])
 a[d..., n]

end

function leading(a::Matrix, n)
 a[1:size(a,1), n]
end

