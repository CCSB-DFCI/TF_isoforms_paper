isos = cloned_isoform_json_array;

function load_isoform_structure(iso) {
    options = {
        customData: { url: 'http://18.210.151.165/tf_isoform/data/alphafold/' + iso['name'] + '/ranked_0.pdb', format: 'pdb' },
        hideControls: "true",
        bgColor: { r: 255, g: 255, b: 255 },
        pdbeLink: false,
        sequencePanel: true
    }
    viewerInstance.visual.update(options, true);
    exon_colors = [];
    for (let i = 0; i < iso.exon_boundaries.length; i++) {
        exon_colors.push({
            struct_asym_id: 'A',
            start_residue_number: iso.exon_boundaries[i][0],
            end_residue_number: iso.exon_boundaries[i][1],
            color: iso.exon_colors[i],
        })
    };
    viewerInstance.events.loadComplete.subscribe(() => {
        viewerInstance.visual.select({
            data: exon_colors,
            nonSelectedColor: { r: 255, g: 255, b: 255 }
        })
    });
};

var viewerInstance = new PDBeMolstarPlugin();
var viewerContainer = document.getElementById('myViewer');
options = {
    customData: { url: 'http://18.210.151.165/tf_isoform/data/alphafold/' + isos[0]['name'] + '/ranked_0.pdb', format: 'pdb' },
    hideControls: "true",
    bgColor: { r: 255, g: 255, b: 255 },
    pdbeLink: false,
    sequencePanel: true
}
viewerInstance.render(viewerContainer, options);


window.onload = function () {
    load_isoform_structure(isos[0]);

    var docFrag = document.createDocumentFragment();
    var makeHandler = function (num) {
        return function () {
            load_isoform_structure(isos[num]);
        };
    };
    for (var i = 0; i < isos.length; i++) {
        var button = document.createElement('button');
        button.innerHTML = isos[i]['name'];
        button.addEventListener('click', makeHandler(i));
        docFrag.append(button);
    }
    document.getElementById('structureButtons').appendChild(docFrag);
};