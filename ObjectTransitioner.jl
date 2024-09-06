function save(name, obj)
    open("object_transfer/$name.txt", "w") do io
        println(io, string(obj))
    end
end

function find(fname)
    open("object_transfer/$fname.txt") do io
        str = read(io, String)
        return parse.(Float64, split(chop(str; head=1, tail=2), ','))
    end
end