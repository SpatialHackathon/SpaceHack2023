# SpaceHack - contributing modules (data, methods, metrics)

Our workflow is set up to allow everyone to contribute "modules" in their preferred programming language (.. as long as that is either R or Python). A module can either be a dataset, a computational method, or an evaluation metric.
![image](https://github.com/SpatialHackathon/SpaceHack2023/assets/114547/7c002916-0a90-4fe7-8745-489313bc0192)

This repository contains some templates and examples of how to implement your module so that it interfaces seamlessly with other modules in the workflow. For example, if you want to implement a new method, you do not need to worry about input data or evaluation metrics as long as you follow the template for reading input and writing output - if you correctly adhere to the input and output guidelines, you should be able to interface with our default data modules and default evaluation metrics modules. The default modules are:
 - data: LIBD Visium DLPFC dataset (4 samples, each with 3 replicates)
 - methods: BayesSpace and SpaGCN
 - evaluation metrics: ARI and V

## Contributing and Code of Conduct

Read our [Contributing Guide](CONTRIBUTING.md) and [Code of Conduct](CODE_OF_CONDUCT.md).

## Contributors

<!-- readme: contributors -start -->
<table>
	<tbody>
		<tr>
            <td align="center">
                <a href="https://github.com/Jieran-S">
                    <img src="https://avatars.githubusercontent.com/u/91852421?v=4" width="75;" alt="Jieran-S"/>
                    <br />
                    <sub><b>Jieran S.</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/shdam">
                    <img src="https://avatars.githubusercontent.com/u/49019552?v=4" width="75;" alt="shdam"/>
                    <br />
                    <sub><b>Søren Helweg Dam</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/niklasmueboe">
                    <img src="https://avatars.githubusercontent.com/u/42138117?v=4" width="75;" alt="niklasmueboe"/>
                    <br />
                    <sub><b>niklasmueboe</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/peicai">
                    <img src="https://avatars.githubusercontent.com/u/55488976?v=4" width="75;" alt="peicai"/>
                    <br />
                    <sub><b>peicai</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/kbiharie">
                    <img src="https://avatars.githubusercontent.com/u/33690856?v=4" width="75;" alt="kbiharie"/>
                    <br />
                    <sub><b>kbiharie</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/heylf">
                    <img src="https://avatars.githubusercontent.com/u/8162688?v=4" width="75;" alt="heylf"/>
                    <br />
                    <sub><b>heylf</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/pakiessling">
                    <img src="https://avatars.githubusercontent.com/u/104848590?v=4" width="75;" alt="pakiessling"/>
                    <br />
                    <sub><b>pakiessling</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Qirongmao97">
                    <img src="https://avatars.githubusercontent.com/u/57286623?v=4" width="75;" alt="Qirongmao97"/>
                    <br />
                    <sub><b>Qirong Mao</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/naveedishaque">
                    <img src="https://avatars.githubusercontent.com/u/114547?v=4" width="75;" alt="naveedishaque"/>
                    <br />
                    <sub><b>Nav</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/zsfrbkv">
                    <img src="https://avatars.githubusercontent.com/u/43470646?v=4" width="75;" alt="zsfrbkv"/>
                    <br />
                    <sub><b>zaira seferbekova</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/markrobinsonuzh">
                    <img src="https://avatars.githubusercontent.com/u/6471769?v=4" width="75;" alt="markrobinsonuzh"/>
                    <br />
                    <sub><b>Mark Robinson</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/sebastiantiesmeyer">
                    <img src="https://avatars.githubusercontent.com/u/25506428?v=4" width="75;" alt="sebastiantiesmeyer"/>
                    <br />
                    <sub><b>sebastiantiesmeyer</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/theinvisibleliya">
                    <img src="https://avatars.githubusercontent.com/u/79532622?v=4" width="75;" alt="theinvisibleliya"/>
                    <br />
                    <sub><b>Liya Zaygerman</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/meghanaturner">
                    <img src="https://avatars.githubusercontent.com/u/22036504?v=4" width="75;" alt="meghanaturner"/>
                    <br />
                    <sub><b>Meghan Turner</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Aokht17">
                    <img src="https://avatars.githubusercontent.com/u/56379827?v=4" width="75;" alt="Aokht17"/>
                    <br />
                    <sub><b>Anastasiia Okhtienko</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/berl">
                    <img src="https://avatars.githubusercontent.com/u/6773896?v=4" width="75;" alt="berl"/>
                    <br />
                    <sub><b>Brian Long</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/gmoranzoni">
                    <img src="https://avatars.githubusercontent.com/u/59561270?v=4" width="75;" alt="gmoranzoni"/>
                    <br />
                    <sub><b>Giorgia Moranzoni</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/tmchartrand">
                    <img src="https://avatars.githubusercontent.com/u/12821536?v=4" width="75;" alt="tmchartrand"/>
                    <br />
                    <sub><b>Tom Chartrand</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/alam-shahul">
                    <img src="https://avatars.githubusercontent.com/u/22669932?v=4" width="75;" alt="alam-shahul"/>
                    <br />
                    <sub><b>alam-shahul</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/svedziok">
                    <img src="https://avatars.githubusercontent.com/u/17719296?v=4" width="75;" alt="svedziok"/>
                    <br />
                    <sub><b>Sven Twardziok</b></sub>
                </a>
            </td>
		</tr>
	<tbody>
</table>
<!-- readme: contributors -end -->

## License

We have adopted the "MIT No Attribution" (MIT-0) License. It is currently attributed to the "SpaceHack organizers", but please also make sure to add your name to your contributions. More on MIT-0 [here](https://github.com/aws/mit-0)
