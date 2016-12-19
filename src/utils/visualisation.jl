#=
    Helper functions to visualise aspects of the algebra.

    Colours were generated using http://colorbrewer2.org/
    --> sequential, 6 classes and then selecting the darkest 3/4

    Plotting is being done using Gadfly:
        http://gadflyjl.org/stable/index.html
=#
const pos_neg_colours = [
    Colors.RGB(a/255,b/255,c/255)
    for (a,b,c) in [
        # blues
        (254,217,118), (219,199,143),
        (49,130,189),(74,119,143),
        (107,174,214),(111,160,176),
        (158,202,225),(143,181,191),
        # purples
        (165,15,21),(135,46,46),
        (117,107,177),(118,113,143),
        (158,154,200),(156,142,173),
        (188,189,220),(180,155,191),
        # reds
        (84,39,143),(76,72,97),
        (222,45,38),(189,68,66),
        (251,106,74),(214,104,71),
        (252,146,114),(217,142,102),
        # greens
        (254,178,76),(219,170,86),
        (49,163,84),(73,140,107),
        (116,196,118),(132,179,148),
        (161,217,155),(159,196,160)
    ]
]

const MUTED_AR_colors = [
    Colors.RGB(a/255,b/255,c/255)
    for (a,b,c) in [
        # blues
        (219,199,143),
        (74,119,143),(111,160,176),(143,181,191),
        # purples
        (135,46,46),
        (118,113,143),(156,142,173),(180,155,191),
        # reds
        (76,72,97),
        (189,68,66),(214,104,71),(217,142,102),
        # greens
        (219,170,86),
        (73,140,107),(132,179,148),(159,196,160)
    ]
]

const AR_colors = [
    Colors.RGB(a/255,b/255,c/255)
    for (a,b,c) in [
        # blues
        (254,217,118),
        (49,130,189),(107,174,214),(158,202,225),
        # purples
        (165,15,21),
        (117,107,177),(158,154,200),(188,189,220),
        # reds
        (84,39,143),
        (222,45,38),(251,106,74),(252,146,114),
        # greens
        (254,178,76),
        (49,163,84),(116,196,118),(161,217,155),
    ]
]

const PALETTE = Scale.color_discrete_manual(AR_colors...)
const MUTED_PALETTE = Scale.color_discrete_manual(MUTED_AR_colors...)


################################
# .: Viewing a Cayley Table :. #
################################
"""
__convert_cayley__

Change the α values in CAYLEY into alternate formats for visualisation
"""
function convert_cayley(output="indices", op=:*, i0=true)
    if i0
        allowed = vcat(ALLOWED[1:end-3], ["10", "20", "30"])
    else
        allowed = vcat(ALLOWED[1:end-3], ["01", "02", "03"])
    end
    cayley = permutedims(
        eval(parse("[α(i) $op α(j) for i in $allowed, j in $allowed]")),
        [1,2]
    )

    if output == "indices"
        return cayley
    elseif output == "strindices"
        return [[a.index for a in cayley[i,:]] for i in 1:16]
    elseif output == "colmap"
        return [[ix(a) for a in cayley[i,:]] for i in 1:16]
    elseif output == "sign"
        return [[a.sign for a in cayley[i,:]] for i in 1:16]
    else
        error("Invalid output specification")
    end
end

"""
__print_cayley__

This is a quick and dirty way to print out the Cayley table in
a couple of different ways so that it can be pasted into excel
for visualising and looking at properties such as sign and symmetry.
"""
function print_cayley(output="indices", op=:*, i0=true, headers=true)
    data = convert_cayley(output, op, i0)
    # TODO:: fix headers for 0i
    headers && println(",$(join(CAYLEY[1,:], ","))")
    for i in 1:16
        c = data[i,:]
        headers && print("$(CAYLEY[1,i]),")
        println(join([a for a in c]..., ","))
    end
end

"""
__visualise_cayley__

Plot a heatmap of the generated Cayley table. If coloured=false then this
will show the disribution of positive and negative products in the algebra.

If a filename is supplied this will save an svg into the current directory.
_NOTE_:: You should not provide a file extension.
"""
function visualise_cayley(;muted=true, coloured=true, dims=(2000,2000),
                          title="Cayley Table for the Williamson Algebra",
                          filename="")
    output = coloured ? "colmap" : "sign"
    pal = muted ? MUTED_PALETTE : PALETTE
    # This is the best way I've found to convert an array of n, n-element
    # arrays into a single nxn multi-array.
    data = hcat(convert_cayley(output)...)

    AR_plot_style = style(
        background_color=colorant"white",
        minor_label_color=colorant"white",
        minor_label_font_size=1px,
        key_position=:none
    )
    Gadfly.push_theme(AR_plot_style)
    plt = spy(
        data,
        Guide.title(title),
        Guide.xlabel(nothing),
        Guide.ylabel(nothing),
        pal
        )

    if filename != ""
        x, y = dims
        if endswith(filename, ".svg")
            img = SVG(filename, x*px, y*px)
        elseif endswith(filename, ".png")
            img = PNG(filename, x*px, y*px)
        elseif endswith(filename, ".pdf")
            img = PDF(filename, x*px, y*px)
        else
            error("filename must end in .svg|png|pdf")
        end
        draw(img, plt)
    else
        return plt
    end
end
