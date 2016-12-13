#=
    Helper functions to visualise aspects of the algebra.

    Colours were generated using http://colorbrewer2.org/
    --> sequential, 6 classes and then selecting the darkest 3/4

    Plotting is being done using Gadfly:
        http://gadflyjl.org/stable/index.html
=#
#=
const MUTED_AR_colors = [
    Colors.RGB(a/255,b/255,c/255)
    for (a,b,c) in [
        # yellow 1
        (219,199,143),
        # reds
        (135,46,46),(189,68,66),(214,104,71),(217,142,102),
        # greens
        (73,140,107),(132,179,148),(159,196,160),
        # blues
        (74,119,143),(111,160,176),(143,181,191),
        # purples
        (118,113,143),(156,142,173),(180,155,191),(76,72,97),
        # yellow 2
        (219,170,86)
    ]
]
=#
const MUTED_AR_colors = [
    Colors.RGB(a/255,b/255,c/255)
    for (a,b,c) in [
        # yellow 1
        (219,199,143), (219,170,86),(76,72,97),(135,46,46),
        # blues
        (74,119,143),(111,160,176),(143,181,191),
        # purples
        (118,113,143),(156,142,173),(180,155,191),
        # reds
        (189,68,66),(214,104,71),(217,142,102),
        # greens
        (73,140,107),(132,179,148),(159,196,160)
    ]
]

const AR_colors = [
    Colors.RGB(a/255,b/255,c/255)
    for (a,b,c) in [
        # yellow 1
        (254,217,118),
        # reds
        (165,15,21),(222,45,38),(251,106,74),(252,146,114),
        # greens
        (49,163,84),(116,196,118),(161,217,155),
        # blues
        (49,130,189),(107,174,214),(158,202,225),
        # purples
        (117,107,177),(158,154,200),(188,189,220),(84,39,143),
        # yellow 2
        (254,178,76)
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
function convert_cayley(output="indices")
    if output == "indices"
        return CAYLEY
    elseif output == "strindices"
        return [[a.index for a in CAYLEY[i,:]] for i in 1:16]
    elseif output == "colmap"
        return [[ix(a) for a in CAYLEY[i,:]] for i in 1:16]
    elseif output == "sign"
        return [[a.sign for a in CAYLEY[i,:]] for i in 1:16]
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
function print_cayley(output="indices", headers=true)
    data = convert_cayley(output)
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
