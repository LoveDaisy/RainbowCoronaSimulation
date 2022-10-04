# Rainbow, Corona Simulation

My simulation program for rainbow and corona.

As we know, tiny water drops in the air may produce rainbow and corona.
Descartes first explain the rainbow with geometry optics over three and a half centuries ago.
Then Airy improved the theory by taking difraction into account one and a half centruies ago.
And half a century later, Lorenz, Mie and Debye derived the most accurate model, which is often called Lorenz-Mie theory and Deby series.

My simulation is based on Lorenz-Mie theory. See [Mie scattering](https://en.wikipedia.org/wiki/Mie_scattering) for detail.
And there is a [draft](Computation%20of%20Mie%20Theory.md) containing all mathematics.

## Rainbow

Size of water drops change the appearance of rainbow dramatically.
When drop size is large, it follows the geometry optics limit, and no supernumerary can be seen;
when drop size is small, diffraction and interference are getting significant, and supernumerary is getting visible.

![rainbow and drop size](img/rainbow_radius_angle.png)

Lee [2] proposed a diagram that described the color changing with water drop size, called Lee diagram.

![rainbow Lee diagram](img/lee_137-145.png)

The colors are rendered with my [color science tools](https://github.com/LoveDaisy/ColorScienceUtils), which converts spectra into RGB colors.

## Secondary rainbow

Colors of secondary rainbow are very similar to primary rainbow, just line up in a reverse order.

![secondary rainbow and drop size](img/sencondary_rainbow_radius_angle.png)

In order to show colors clearly, I normalize them by max illuminance. They are in fact much fainter than primary rainbow.

Also here is the Lee diagram for secondary rainbow,

![secondary rainbow Lee diagram](img/lee_123-131.png)

## Corona

Corona is often produced by high clouds. The particles in high clouds are much smaller than rain drops. This leads to more diffraction effects than case of rainbow.
Similarly, with particle size goes smaller, the colors become more spread.

![corona and drop size](img/corona_radius_angle.png)

To display the faint colors of corona clearly, I amplify the intensity by factor of 4, so the bright parts are over saturated to white. The same to the Lee diagram.

![corona Lee diagram](img/lee_000-010.png)

Suprisingly, there are *stripes* on Lee diagram when particle size is very small.

## Reference

[1] Laven, Philip. "Simulation of rainbows, coronas, and glories by use of Mie theory." Applied optics 42.3 (2003): 436-444.

[2] Lee, Raymond L. "Mie theory, Airy theory, and the natural rainbow." Applied Optics 37.9 (1998): 1506-1519.

[3] Wang, Ru T., and H. C. Van de Hulst. "Rainbows: Mie computations and the Airy approximation." Applied optics 30.1 (1991): 106-117.

[4] [Philip Laven's website](http://www.philiplaven.com/index1.html) and his [MiePlot](http://www.philiplaven.com/mieplot.htm) program.
