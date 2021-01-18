
document.addEventListener('toggled', function(event) {
  var status = document.getElementById('shiny-toggle-status');
  if (event.now_active === 'before') {
    status.innerHTML = 'Before';
  } else {
    status.innerHTML = 'After';
  }
});

// Disable/Enable validate buttons
Shiny.addCustomMessageHandler('toggle-validate-btns-handler', toggleBtns);
function toggleBtns(message) {
  document.getElementById('group_validation_button').disabled = message;
  document.getElementById('case_validation_button').disabled = message;
}

//Keyboard shortcuts
document.onkeydown = function(e) {
  if (e.shiftKey) { // using e.shiftKey && e.which == 86 simultaneously doesn't work?
    if (e.which == 13) { // ENTER
      Shiny.onInputChange("validateGroup", Math.random());
    }
  } else {
     if (e.which == 39) { // right-arrow
      Shiny.onInputChange("nextCase", Math.random());
     } else if (e.which == 37) { // left-arrow
      Shiny.onInputChange("prevCase", Math.random());
     } else if (e.which == 40) { // down-arrow
      Shiny.onInputChange("nextType", Math.random());
     } else if (e.which == 38) { //up-arrow
      Shiny.onInputChange("prevType", Math.random());
     } else if (e.which == 13) { // ENTER
      Shiny.onInputChange("validateCase", Math.random());
    } else if (e.which == 27) { //esc, q = 81
      Shiny.onInputChange("quit_button", Math.random());
    }
  }
};
