r = Corpus(ReutersSource("reut2-000.xml"))
commodity = c("gold","copper","tin","alum","iron-steel","lead","nickel","palladium",
  "platinum","silver","strategic-metal","tung","zin","carcass","livestock","meal-feed",
  "rubber","barley","cocoa","coffee","cotton","grain","oilseed","soybean","tea","sugar","wheat","veg-oil")
energy = c("crude","nat-gas","heat","propane","fuel","gas","pet-chem","jet","naptha")
econ = c("bop","cpi","jobs","gnp","reserves","ipi","lei","housing","trade","money-supply","retail","money-fx","interest","ship")
currency = c("dlr","austdlr","hk","singdlr","nzdlr","can","stg","dmk","yen","sfr")
corp=c("acq","earn")
money.fx = c("money-fx")
ship = c("ship")
interest = c("interest")

#misc = c("money-fx","ship","interest")

#x.topics[x.topics %in% metals] = 1 # black 
x.topics[x.topics %in% commodity] = "black" #
x.topics[x.topics %in% energy] = "red" #red
x.topics[x.topics %in% econ] = "green" #green
x.topics[x.topics %in% corp] = "blue" #blue
x.topics[x.topics %in% money.fx] = "yellow" #pink
x.topics[x.topics %in% ship] = "pink" #pink
x.topics[x.topics %in% interest] = "purple"
#N = nrow(xyz)

N = 1000
topics = seq("N/A", "N/A", length.out = N)
# Get the topics list
for(i in 1:N){
  ri.meta = meta(r[[i]])
  xyz = ri.meta$Topics
  if(length(xyz) > 0){
    if(length(xyz) == 1){
      if(xyz %in% commodity){
        topics[i] = "black"
      }
      else if(xyz %in% energy){
        topics[i] = "red" 
      }
      else if(xyz %in% econ){
        topics[i] = "green" 
      }
      else if(xyz %in% corp){
        topics[i] = "blue" 
      }
      else if(xyz %in% money.fx){
        topics[i] = "yellow"
      }
      else if(xyz %in% ship){
        topics[i] = "pink"
      }
      else if(xyz %in% interest){
        topics[i] = "purple"
      }
      else{
        topics[i] = "magenta" 
      }
    }
    if(length(xyz) > 1){
      xyz[xyz %in% commodity] = "black" 
      xyz[xyz %in% energy] = "red" 
      xyz[xyz %in% econ] = "green" 
      xyz[xyz %in% corp] = "blue" 
      xyz[xyz %in% money.fx] = "yellow" 
      xyz[xyz %in% ship] = "pink" 
      xyz[xyz %in% interest] = "purple" 
      freqs = table(xyz)
      ind = order(freqs,decreasing = TRUE)
      topics[i] = names(freqs)[ind[1]]
      if(length(freqs >= 2)){
        topics[i] = "magenta"
       }
    }
  }
}

index = c(1:N)
index = index[topics != "N/A"]
r.new = r[index]

reuters = tm_map(r.new, as.PlainTextDocument)
reuters = tm_map(reuters, stripWhitespace)
reuters = tm_map(reuters,tolower)
reuters = tm_map(reuters,removeWords,stopwords("english"))
reuters = tm_map(reuters, stemDocument)
dtm = DocumentTermMatrix(reuters)
dtm.sparse = removeSparseTerms(dtm, 0.95)
dtm.mat = as.matrix(dtm.sparse)
W = cosine(t(dtm.mat))
P = transition.matrix(W)
## Embed using expected commute time
dist.ect = ect(P)
loc = cmdscale(dist.ect, k = 2)
plot(loc,col=x.topics)


smartlegend(x="left",y="top", inset = 0,c("commodity", "energy","economic","corporation","miscellaneous"),fill = c("black","red","green","blue","pink"))

rng_sphere <- function(d, type='rgl')
{

n <- length(d)
nn <- n - 3

# init our x,y,z coordinate arrays
x <- array(dim=nn)
y <- array(dim=nn)
z <- array(dim=nn)

# init red,green,blue color component arrays
cr <- array(dim=nn)
cg <- array(dim=nn)
cb <- array(dim=nn)


# convert lagged random numbers from d into spherical coordinates
# then convert to cartesian x,y,z coordinates for simple display
for (i in 1:nn)
{
theta <- 2*pi*d[i]
phi <- pi*d[i+1]
r <- sqrt(d[i+2])

x[i] <- r * sin(theta) * cos(phi)
y[i] <- r * sin(theta) * sin(phi)
z[i] <- r * cos(theta)

# give each location a color based on some rules
cr[i] <- d[i] / max(d)
cg[i] <- d[i+1] / max(d)
cb[i] <- d[i+2] / max(d)

} # end function


if( type == 'rgl')
{
# setup rgl environment:
zscale <- 1
 
# clear scene:
clear3d("all")
 
# setup env:
bg3d(color="white")
light3d()
 
# draw shperes in an rgl window
spheres3d(x, y, z, radius=0.025, color=rgb(cr,cg,cb))
}

if(type == '2d')
{
# optional scatterplot in 2D
scatterplot3d(x,y,z, pch=16, cex.symbols=0.25, color=rgb(cr,cg,cb), axis=FALSE )
}


# optionally return results
return (list(x=x, y=y, z=z, red=cr, green=cg, blue=cb))

}



