var idChecked = false;


jQuery(document).ready(function() {

jQuery('a.login-window').click(function() {
    
    //Getting the variable's value from a link 
    var loginBox = jQuery(this).attr('href');

    //Fade in the Popup
    jQuery(loginBox).fadeIn(300);
    
    //Set the center alignment padding + border see css style
    var popMargTop = (jQuery(loginBox).height() + 24) / 2; 
    var popMargLeft = (jQuery(loginBox).width() + 24) / 2; 
    
    jQuery(loginBox).css({ 
        'margin-top' : -popMargTop,
        'margin-left' : -popMargLeft
    });
    
    // Add the mask to body
    jQuery('body').append('<div id="mask"></div>');
    jQuery('#mask').fadeIn(300);
    
    return false;
});

jQuery(window).keydown(function(event) {
  if(!idChecked)
    if (event.keyCode == 13) {
      event.preventDefault();
      //jQuery("form").submit();
      checkIds2();
      return false;
    }
});


});
