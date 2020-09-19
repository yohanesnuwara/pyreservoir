# PyReservoir
Python utilities for reservoir engineering calculations

<div>
<img src="https://user-images.githubusercontent.com/51282928/85827088-bb6f1300-b7af-11ea-9a1f-eed08adddaff.png" width="500"/>
</div>

* PVT analysis `pvt`
* Volumetric mapping `volumetrics`
* Well-test analysis `welltest`
* Material balance `matbal`
* Decline curve analysis `dca`

## Applications

* `pvt`: Obtaining PVT properties using PVT correlation ([`pvtcorrelation`](https://github.com/yohanesnuwara/pyreservoir/blob/master/pvt/pvtcorrelation.py)); Processing PVT lab experiments such as DL, CCE, Flash Analysis, and CVD ([`pvtlab`](https://github.com/yohanesnuwara/pyreservoir/blob/master/pvt/pvtlab.py))
* `volumetrics`: Computing OOIP and OGIP using volumetric trapezoidal, pyramidal, and Simpson's 1/3 rule ([`volumetrics`](https://github.com/yohanesnuwara/pyreservoir/blob/master/volumetrics/volumetrics.py))
* `welltest`
* `matbal`: Computing aquifer influx into a reservoir using Schilthuis, van Everdingen-Hurst, and Fetkovich methods ([`aquifer`](https://github.com/yohanesnuwara/pyreservoir/blob/master/matbal/aquifer.py)); Material balance plots for OOIP and OGIP verification for gas and oil reservoirs ([`mbal`](https://github.com/yohanesnuwara/pyreservoir/blob/master/matbal/mbal.py)); Computing reservoir drive indices and producing energy plots ([`drives`](https://github.com/yohanesnuwara/pyreservoir/blob/master/matbal/drives.py))
* `dca`: Computing decline-curve parameters from Arps 

## Tutorials

See these tutorial [Jupyter notebooks](https://github.com/yohanesnuwara/pyreservoir/tree/master/notebooks) to get started with each of the libraries. See inside [this folder](https://github.com/yohanesnuwara/pyreservoir/tree/master/data) to get the datasets used for the tutorials. 

## License

I consider the goodness of open-source program but I strongly recommend that anyone who wish to use any program in this package to consider the code authorship. This work is licensed with Creative Commons BY-NC-ND 4.0 International. 

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.
