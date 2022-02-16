#containers for holding information about our analysis

struct member

    node_i::Int64
    node_j::Int64
    cross_section::String
    material::String

end

struct node

    x::Float64
    y::Float64

end


struct material

    name::String
    E::Float64

end

struct cross_section

    name::String
    A::Float64

end
##
#define those containers
#user defined...

materials = [material("steel", 29000.0), material("aluminium", 15000.0), material("concrete", 4000)]

nodes = [node(0.0, 0.0), node(9.0*12, 11.0*12), node(18.0*12, 0.0)]

cross_sections = [cross_section("I-beam", 20.0), cross_section("angle", 10.0)]

members = [member(1, 2, "I-beam", "steel"), member(2, 3, "angle", "steel")]

##
# global varaibles
global F1 = 0
global F2 = 0
global x = 0
global y = 0

global distortion = [x
                     y]
global forces = [F1
                 F2]

global overlap = [0 0
                  0 0]
global value = 1

##

#define the local element stiffness matrix
function construct_local_element_k_matrix(E, A, L)
    println()
    print("Local ")
    println(value)


    k = (E*A/L)* [1.0 0.0 -1.0 0.0
                 0.0 0.0 0.0 0.0
               -1.0 0.0 1.0 0.0
               0.0 0.0 0.0 0.0]





    print_symetric_matrix(round_matrix(k))
    println()

    global value = value + 1

    return k

end
## Rounding Function
function round_matrix(matrix)
    a = (size(matrix)[1])^2
    b = size(matrix)[1]
    i = 1

    rounded_matrix = zeros(b, b)

    while  i <= a
        rounded_matrix[i] = round(matrix[i], digits=2)
        i = i + 1
    end
    return rounded_matrix
end

##Printing function

function print_symetric_matrix(matrix)
  rownumber = 1
  a = size(matrix)[1]
  b = 1
  c = 0
  i = 1

  while rownumber <= a
    print("[")

    while i <= a
      print(matrix[b])

      if i < a
        print(", ")
        b = b + a
      else
        b = b + a
        println("]")
      end

      i = i + 1
    end

    i = i - a
    rownumber = rownumber + 1
    b = rownumber
  end
end



##
#Calculations

    # find member length
function find_member_length(node_i, node_j)
        x1 = node_i.x
        x2 = node_j.x
        y1 = node_i.y
        y2 = node_j.y

        member_length = sqrt((x1 - x2)^2 + (y1 - y2)^2)

    return member_length
end
##Matrix Transformations

function construct_xy_rotation_matrix(angle, matrix)
    angle = deg2rad(angle)
        T_matrix =       [cos(angle) -sin(angle) 0 0
                          sin(angle) cos(angle) 0 0
                          0 0 1 0
                          0 0 0 1]
        xy_matrix =  (T_matrix * matrix) * inv(T_matrix)
    return xy_matrix
end

function construct_yz_rotation_matrix(angle, matrix)
    angle = deg2rad(angle)
            T_matrix = [1 0 0 0
                        0 cos(angle) -sin(angle) 0
                        0 sin(angle) cos(angle) 0
                        0 0 0 1]
        yz_matrix =  (T_matrix * matrix) * inv(T_matrix)
        return yz_matrix
end

function construct_xz_rotation_matrix(angle, matrix)
    angle = deg2rad(angle)
              T_matrix = [cos(angle) 0 sin(angle) 0
                          0 1 0 0
                          -sin(angle) 0 cos(angle) 0
                          0 0 0 1]
        xz_matrix =  (T_matrix * matrix) * inv(T_matrix)
        return xz_matrix
end

#I think something is wrong with this one but im not sure what
function construct_other_full_rotation_matrix(ang, matrix)
    angle = deg2rad(ang)
    T_matrix =           [cos(angle) -sin(angle) 0 0
                          sin(angle) cos(angle) 0 0
                          0 0 cos(angle) -sin(angle)
                          0 0 sin(angle) cos(angle)]

        full_matrix =  T_matrix * matrix * inv(T_matrix)


        print("T_matrix ")
        println(value -1)
    print_symetric_matrix(round_matrix(full_matrix))
    println()
        return full_matrix
end

function construct_full_rotation_matrix(ang, matrix)
    angle = deg2rad(ang)
    cos2 = (cos(angle))^2
    cosSin = cos(angle) * sin(angle)
    sin2 = (sin(angle))^2

    T_matrix =  matrix * [cos2 cosSin -cos2 -cosSin
                        cosSin sin2 -cosSin -sin2
                        -cos2 -cosSin cos2 cosSin
                        -cosSin -sin2 cosSin sin2]

    print("T_matrix ")
    println(value -1)
    print_symetric_matrix(round_matrix(T_matrix))
    println()

return T_matrix

end

function print_invert_matrix(matrix)
    inverted_matrix = inv(matrix)
    println("Inverted")
    print_symetric_matrix(round_matrix(inverted_matrix))

end

##  global construction

function construct_global_matrix_2_elements(m1, m2)

    m_global = [m1[1] m1[5] m1[9] m1[13] 0 0
               m1[2] m1[6] m1[10] m1[14] 0 0
               m1[3] m1[7] (m1[11]+m2[1]) (m1[15]+m2[5]) m2[9] m2[13]
               m1[4] m1[8] (m1[12]+m2[2]) (m1[16]+m2[6]) m2[10] m2[14]
               0 0 m2[3] m2[7] m2[11] m2[15]
               0 0 m2[4] m2[8] m2[12] m2[16]]
    println("Global")
    print_symetric_matrix(round_matrix(m_global))
    println()


    m_overlap = [m_global[15] m_global[21]
                 m_global[16] m_global[22]]


    println("Overlap")
    print_symetric_matrix(round_matrix(m_overlap))
    println()
    global overlap= m_overlap
    return (m_overlap, m_global)
end



## Commands

construct_global_matrix_2_elements(
    construct_full_rotation_matrix(
        41,
        construct_local_element_k_matrix(
            materials[1].E,
            cross_sections[1].A,
            find_member_length(nodes[1], nodes[2])
        )
    ),

    construct_full_rotation_matrix(
        48,
        construct_local_element_k_matrix(
            materials[1].E,
            cross_sections[2].A,
            find_member_length(nodes[2], nodes[3])
        )
    )
)


