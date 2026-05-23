# Install htmlwidgets if you haven't already:
# install.packages("htmlwidgets")

library(networkD3)
library(htmlwidgets) # This lets us add custom JavaScript

# 1. Define the Nodes
nodes <- data.frame(
  id = 0:5, 
  name = c("Placebo", "Adalimumab", "Ustekinumab", "Secukinumab", "Ixekizumab", "Guselkumab"),
  group = c("Placebo", "TNF-alpha", "IL-12/23", "IL-17", "IL-17", "IL-23"),
  participants = c(2500, 1800, 1200, 1500, 900, 800)
)
nodes$node_size <- sqrt(nodes$participants) * 1.5

# 2. Define the Links
links <- data.frame(
  source = c(0,  0,  0,  0,  0,  1,  2,  3),
  target = c(1,  2,  3,  4,  5,  3,  4,  5),
  num_studies = c(15, 8, 12, 6,  4,  3,  2,  1) 
)

# 3. Generate the Base Network Diagram
nma_plot <- forceNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "num_studies",   
  NodeID = "name",         
  Nodesize = "node_size",  
  Group = "group",         
  opacity = 0.9,  
  opacityNoHover = 1,
  zoom = TRUE,             
  fontSize = 14,
  fontFamily = "Arial",
  linkDistance = 150,      
  charge = -400,           # Slightly increased negative charge to make room for text
  bounded = TRUE,          
  legend = TRUE            
)

# 4. Inject JavaScript to force permanent labels and add non-destructive hover effects
nma_plot_labelled <- onRender(nma_plot, '
  function(el, x) {
    // 1. Remove the default networkD3 hover text and titles
    d3.selectAll(".node text").remove();
    d3.selectAll(".node title").remove(); 
    
    // 2. Append our permanent labels next to the nodes
    d3.selectAll(".node").append("text")
      .attr("dx", 20)         
      .attr("dy", ".35em")    
      .style("font-size", "14px")
      .style("font-family", "Arial")
      .style("transition", "all 0.2s ease-in-out") 
      .style("pointer-events", "none") 
      .text(function(d) { return d.name; });
      
    // 3. Add namespaced mouse events (.labelHover) so we don\'t break default networkD3 highlighting
    d3.selectAll(".node")
      .on("mouseover.labelHover", function(d) {
        d3.select(this).select("text")
          .style("font-size", "18px")
          .style("font-weight", "bold");
      })
      .on("mouseout.labelHover", function(d) {
        d3.select(this).select("text")
          .style("font-size", "14px")
          .style("font-weight", "normal");
      });
  }
')

# View the perfectly labeled plot
nma_plot_labelled


saveWidget(nma_plot_labelled, file = "C:/Users/p12916tc/Documents/GitHub/nma/diagram.html")
