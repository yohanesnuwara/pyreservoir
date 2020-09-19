# PyReservoir
Python utilities for reservoir engineering calculations

<div>
<img src="https://user-images.githubusercontent.com/51282928/85827088-bb6f1300-b7af-11ea-9a1f-eed08adddaff.png" width="500"/>
</div>

* PVT analysis `pvt`
* Volumetric mapping `pyvolumetric`
* Well-test analysis `welltest`
* Material balance `matbal`
* Decline curve analysis `dca`

## MBAL for dry-gas (DG) and gas-condensate (GC)

**Differences and similarities of variables**

|**Variables**|**DG** & **GC above dew-point**|**GC below dew-point**|
|:--:|:--:|:--:|
|Reservoir voidage (<img src="https://render.githubusercontent.com/render/math?math=F">)|<img src="https://render.githubusercontent.com/render/math?math=F=G_pB_g">|<img src="https://render.githubusercontent.com/render/math?math=F=[N_p(\frac{B_o-R_sB_g}{1-R_vR_s})]+[(G_p-G_i)(\frac{B_g-R_vB_o}{1-R_vR_s})]">|
|Total gas FVF (<img src="https://render.githubusercontent.com/render/math?math=B_{tg}">)|<img src="https://render.githubusercontent.com/render/math?math=B_{tg}=B_g">|<img src="https://render.githubusercontent.com/render/math?math=B_{tg}=\frac{B_g(1-R_vR_{vi})+(R_{vi}-R_v)B_o}{1-R_vR_s}">|
|Gas expansion factor (<img src="https://render.githubusercontent.com/render/math?math=E_g">)|<img src="https://render.githubusercontent.com/render/math?math=E_g=B_{tg}-B_{gi}">|<img src="https://render.githubusercontent.com/render/math?math=E_g=B_{tg}-B_{gi}">|

## License

I consider the goodness of open-source program but I strongly recommend that anyone who wish to use any program in this package to consider the code authorship. This work is licensed with Creative Commons BY-NC-ND 4.0 International. 

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.
