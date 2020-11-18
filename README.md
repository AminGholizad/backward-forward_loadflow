# backward-forward_loadflow
Backward/forward loadflow of radial distribution systems.
Inputs are line data, node number of the feader start (root node), base voltage, and load data of the system, any additional destributed generations or reactive compensators.
Loads are considered as constant power.

[![View backward-forward_loadflow on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/70024-backward-forward_loadflow)

# Examples
There are two examples of ieee33bus and ieee69bus systems that will plot voltage profile of the loadflow results.

# TODO
Add support for other load types (constant current, constant impedance)
