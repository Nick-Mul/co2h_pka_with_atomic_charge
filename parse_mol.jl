using JSON
using LinearAlgebra

function read_json(filename)
    # Directly parse the JSON file
    data = JSON.parsefile(filename)
    
    return data
end

function read_mol2(filename::String)
    lines = readlines(filename)
    atoms = []
    bonds = []
    charges = []
    reading_atoms = false
    reading_bonds = false

    for line in lines
        if startswith(line, "@<TRIPOS>ATOM")
            reading_atoms = true
            reading_bonds = false
            continue
        elseif startswith(line, "@<TRIPOS>BOND")
            reading_atoms = false
            reading_bonds = true
            continue
        elseif startswith(line, "@<TRIPOS>")
            reading_atoms = false
            reading_bonds = false
            continue
        end

        if reading_atoms
            parts = split(line)
            if length(parts) >= 9
                push!(atoms, (parts[1], parts[2], parse(Float64, parts[3]), parse(Float64, parts[4]), parse(Float64, parts[5]), parts[6], parts[7], parts[8], parse(Float64, parts[9])))
            end
        elseif reading_bonds
            parts = split(line)
            if length(parts) >= 4
                push!(bonds, (parts[1], parts[2], parts[3], parts[4]))
            end
        end
    end

    return atoms, bonds
end

function write_mol2(filename::String, atoms, bonds)
    open(filename, "w") do io
        println(io, "@<TRIPOS>MOLECULE")
        println(io, "Molecule Name")
        println(io, "$(length(atoms)) $(length(bonds))")
        println(io, "SMALL")
        println(io, "USER_CHARGES")
        println(io, "\n")

        println(io, "@<TRIPOS>ATOM")
        for (i, atom) in enumerate(atoms)
            println(io, "$(atom[1])  $(atom[2])  $(atom[3])  $(atom[4])  $(atom[5])  $(atom[6])  $(atom[7])  $(atom[8])  $(atom[9])")
        end

        println(io, "\n@<TRIPOS>BOND")
        for (i, bond) in enumerate(bonds)
            println(io, "$(bond[1])  $(bond[2])  $(bond[3])  $(bond[4])")
        end
    end
    println("MOL2 file written to $filename")
end

function add_charges(atoms, new_charges)
    if length(atoms) != length(new_charges)
        error("Number of atoms and charges must match")
    end
    
    new_atoms = []
    for (i, atom) in enumerate(atoms)
        new_atom = (atom[1], atom[2], atom[3], atom[4], atom[5], atom[6], atom[7], atom[8], new_charges[i])
        push!(new_atoms, new_atom)
    end
    
    return new_atoms
end

file = ARGS[1]
base, filetype = split(file, ".")
print(base)
run(`obabel $base.xyz -omol2 -O $base.mol2`)
atoms, bonds = read_mol2("$base.mol2")
data = read_json("$base.property.json")
nested_list = collect(data["Geometry_1"]["Mulliken_Population_Analysis"]["ATOMICCHARGES"])
new_charges = [item[1] for item in nested_list]
print(new_charges)
print(atoms)
atoms_with_charges = add_charges(atoms, new_charges)
output_file = "$base.mol2"
write_mol2(output_file, atoms_with_charges, bonds)
