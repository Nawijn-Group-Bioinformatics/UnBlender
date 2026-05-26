Shiny.addCustomMessageHandler("testmessage",
  function(message) {
    alert(JSON.stringify(message));
  }
);

function hide_all_abstracts(mode){
    console.log('toggling all')
    var myClasses = document.querySelectorAll('.medline_abstract_dynamic'),
        i = 0,
        l = myClasses.length;

    for (i; i < l; i++) {
        if (mode =='collapse'){
            myClasses[i].style.display = 'none';
        }
        else{
            myClasses[i].style.display = 'block';
        }
    }
    var myClasses = document.querySelectorAll('.abstract_show'),
        i = 0,
        l = myClasses.length;

    for (i; i < l; i++) {
        if (mode =='collapse'){
            myClasses[i].style.display = 'inline';
        }
        else{
            myClasses[i].style.display = 'none';
        }
    }
    var myClasses = document.querySelectorAll('.abstract_hide'),
            i = 0,
            l = myClasses.length;

    for (i; i < l; i++) {
        if (mode =='collapse'){
            myClasses[i].style.display = 'none';
        }
        else{
            myClasses[i].style.display = 'inline';
        }
    }

};





function toggle(id){
   console.log('toggling single')
   div_el  = document.getElementById(id);
   olink_el = document.getElementById('olink_'+id);
   hlink_el = document.getElementById('hlink_'+id);
   if (div_el.style.display != 'none')
   {

      div_el.style.display='none';
      olink_el.style.display='inline';
      hlink_el.style.display='none';
      /*img_el.src = '../show.png';*/


   }
   else
   {
      div_el.style.display='inline';
      olink_el.style.display='none';
      hlink_el.style.display='inline';



   };
};



