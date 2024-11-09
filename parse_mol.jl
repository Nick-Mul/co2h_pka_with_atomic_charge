using JSON
using LinearAlgebra

function read_json(filename)
    open(filename, "r") do io
        return JSON.parse(io)
    end
end

function read_mol2(filename::String)
    atoms = Vector{NamedTuple{(:id, :name, :x, :y, :z, :type, :residue_id, :residue_name, :charge), Tuple{String, String, Float64, Float64, Float64, String, String, String, Float64}}}()
    bonds = Vector{NamedTuple{(:id, :atom1, :atom2, :type), Tuple{String, String, String, String}}}()
    
    open(filename, "r") do io
        section = :none
        for line in eachline(io)
            if startswith(line, "@<TRIPOS>ATOM")
                section = :atom
            elseif startswith(line, "@<TRIPOS>BOND")
                section = :bond
            elseif startswith(line, "@<TRIPOS>")
                section = :none
            elseif section == :atom
                parts = split(line)
                if length(parts) >= 9
                    push!(atoms, (id=parts[1], name=parts[2], x=parse(Float64, parts[3]), y=parse(Float64, parts[4]), 
                                  z=parse(Float64, parts[5]), type=parts[6], residue_id=parts[7], residue_name=parts[8], 
                                  charge=parse(Float64, parts[9])))
                end
            elseif section == :bond
                parts = split(line)
                if length(parts) >= 4
                    push!(bonds, (id=parts[1], atom1=parts[2], atom2=parts[3], type=parts[4]))
                end
            end
        end
    end
    
    return atoms, bonds
end


"""
Write molecular structure data to a MOL2 file format.
"""
function write_mol2(filename::String, atoms, bonds)
    open(filename, "w") do io
        println(io, "@<TRIPOS>MOLECULE")
        println(io, "Molecule Name")
        println(io, "$(length(atoms)) $(length(bonds))")
        println(io, "SMALL")
        println(io, "USER_CHARGES")
        println(io, "\n")

        println(io, "@<TRIPOS>ATOM")
        for atom in atoms
            println(io, "$(atom.id)  $(atom.name)  $(atom.x)  $(atom.y)  $(atom.z)  $(atom.type)  $(atom.residue_id)  $(atom.residue_name)  $(atom.charge)")
        end

        println(io, "\n@<TRIPOS>BOND")
        for bond in bonds
            println(io, "$(bond.id)  $(bond.atom1)  $(bond.atom2)  $(bond.type)")
        end
    end
    println("MOL2 file written to $filename")
end


"""
Add new charges to a list of atoms, replacing their existing charges.
# returns
A new vector of NamedTuples with the same structure as the input `atoms`, but with the `charge` field updated
to the corresponding value from `new_charges`.
"""
function add_charges(atoms, new_charges)
    length(atoms) == length(new_charges) || error("Number of atoms and charges must match")
    
    return [NamedTuple{(:id, :name, :x, :y, :z, :type, :residue_id, :residue_name, :charge)}(
        (atom.id, atom.name, atom.x, atom.y, atom.z, atom.type, atom.residue_id, atom.residue_name, charge)
    ) for (atom, charge) in zip(atoms, new_charges)]
end

function main(args)
    isempty(args) && error("Please provide a file name as an argument")
    file = args[1]
    base, _ = splitext(file)
    
    run(`obabel $base.xyz -omol2 -O $base.mol2`)
    atoms, bonds = read_mol2("$base.mol2")
    data = read_json("$base.property.json")
    new_charges = first.(data["Geometry_1"]["Mulliken_Population_Analysis"]["ATOMICCHARGES"])
    
    atoms_with_charges = add_charges(atoms, new_charges)
    output_file = "$base.mol2"
    write_mol2(output_file, atoms_with_charges, bonds)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end