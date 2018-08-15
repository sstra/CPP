
role <- function(r) {
  span(class = 'label label-primary', r)
}




tabAbout = tabPanel("About",
  div(class = "panel panel-default",
    div(class = 'panel-heading', 
        h3(class='panel-title', style = 'font-weight: bold', 
           "Primary Contributor")
       ),
    div(class = "panel-body",
        p(span(style='font-weight:bold', "Garrett M. Dancik, PhD"),
         "is an Associate Professor of Computer Science / ",
         a(href = "https://gdancik.github.io/bioinformatics/", "Bioinformatics"), "at Eastern Connecticut State University (Willimantic, CT).", 
         role("Maintainer"),
         role("Contributer")
         ) 
    ),

    div(class = 'panel-heading', h3(class='panel-title', style = 'font-weight: bold', "Additional Contributors")),

    div(class = "panel-body",
          p("Person 1",role("Contributor")),
          p("Person 2", role("Contributor")),
          p("Stefanos Stravoravdis is a biology and mathematics major at Eastern Connecticut State University (Willimantic, Connecticut)",
            "with research interests involving bioinformatics and studying resistance in microorganisms.",
            role("Contributor"))
    )
  )
)

