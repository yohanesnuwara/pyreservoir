# PyReservoir
Python utilities for reservoir engineering calculations

<div>
<img src="https://user-images.githubusercontent.com/51282928/85827088-bb6f1300-b7af-11ea-9a1f-eed08adddaff.png" width="500"/>
</div>

* PVT analysis `pypvt`
* Volumetric mapping `pyvolumetric`
* Well-test analysis `pywelltest`
* Material balance `pymbal`
* Decline curve analysis `pydca`

## MBAL for dry-gas (DG) and gas-condensate (GC)

**Differences and similarities of variables**

|**Variables**|**DG** & **GC above dew-point**|**GC below dew-point**|
|:--:|:--:|:--:|
|Reservoir voidage (<img src="https://render.githubusercontent.com/render/math?math=F">)|<img src="https://render.githubusercontent.com/render/math?math=F=G_pB_g">|<img src="https://render.githubusercontent.com/render/math?math=F=N_p(\frac{B_o-R_sB_g}{1-R_vR_s})+(G_p-G_i)(\frac{B_g-R_vB_o}{1-R_vR_s})">|
|Total gas FVF (<img src="https://render.githubusercontent.com/render/math?math=B_{tg}">)|<img src="https://render.githubusercontent.com/render/math?math=B_{tg}=B_g">|<img src="https://render.githubusercontent.com/render/math?math=B_{tg}=\frac{B_g(1-R_vR_{vi})+(R_{vi}-R_v)B_o}{1-R_vR_s}">|
|Gas expansion factor (<img src="https://render.githubusercontent.com/render/math?math=E_g">)|<img src="https://render.githubusercontent.com/render/math?math=E_g=B_{tg}-B_{gi}">|<img src="https://render.githubusercontent.com/render/math?math=E_g=B_{tg}-B_{gi}">|
