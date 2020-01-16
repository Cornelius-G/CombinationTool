function findfirstduplicate(x::AbstractArray{T}) where T
    uniqueset = Set{T}()
    for i in x
        if(i in uniqueset)
            return true, i
        end
        push!(uniqueset, i)
    end
    return false
end


function findduplicates(x::AbstractArray{T}) where T
    uniqueset = Set{T}()
    duplicates = Vector{T}()
    duplicate_idx = Vector{T}()

    for i in x
        if(i in uniqueset)
            push!(duplicates, i)
        else
            push!(uniqueset, i)
        end
    end
    return duplicates, duplicate_idx
end


→(val, s) = getfield(val, s)
