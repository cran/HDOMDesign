\name{HDOMdesign}
\alias{HDOMdesign}
\alias{HDOM2}
\alias{HDOM4}
\alias{HDOM8}
\alias{HDOM16}
\title{High-dimensional orthogonal maximin distance designs}
\description{
Generates a high-dimensional orthogonal maximin distance designs. 
}
\usage{
HDOMdesign(n,p,s)
HDOM2(n,p)
HDOM4(n,p)
HDOM8(n,p)
HDOM16(n,p)
}
\arguments{
\item{n}{Number of points.}
\item{p}{Number of dimensions.}
\item{s}{Number of levels.}
}
\details{
These functions generate high-dimensional orthogonal maximin distance designs. 
In particular, HDOM2, HDOM4, HDOM8, and HDOM16 generate designs in two, four, eight, and sixteen levels.
}
\value{
The value returned from the function is an n times p matrix giving the experimental design}
%\references{}
\examples{
HDOM2(n=80,p=60)
HDOM4(n=80,p=60)
HDOM8(n=80,p=60)
HDOM16(n=80,p=60)
HDOMdesign(n=80,p=60,s=4)
}
