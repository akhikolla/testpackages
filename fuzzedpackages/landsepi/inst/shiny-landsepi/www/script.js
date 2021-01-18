setTimeout(function() {
    /********************************************************************/
    $("#cultivars").on("click", "td", function() {
        var table = $('#cultivars table').DataTable();
        if (table.cell(this).data() === 0) {
            table.cell(this).data(1);
        } else if (table.cell(this).data() === 1) {
            table.cell(this).data(0);
        }
    });
    /********************************************************************/
    $("#croptypes").on("click", "td", function() {
        
    });
    /********************************************************************/
}, 0);

