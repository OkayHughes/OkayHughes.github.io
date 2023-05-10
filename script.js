/*
This is your site JavaScript code - you can add interactivity and carry out processing
- Initially the JS writes a message to the console, and moves a button you can add from the README
*/

// Print a message in the browser's dev tools console each time the page loads
// Use your menus or right-click / control-click and choose "Inspect" > "Console"
console.log("Hello 🌎");

/* 
Make the "Click me!" button move when the visitor clicks it:
- First add the button to the page by following the "Next steps" in the README
*/
// const btn = document.querySelector("button"); // Get the button from the page
// // Detect clicks on the button
// if (btn) {
//   btn.onclick = function() {
//     // The JS works in conjunction with the 'dipped' code in style.css
//     btn.classList.toggle("dipped");
//   };
// }


(function() {
  const width = window.innerWidth;
  const height = window.innerHeight;
  const ratio = (width / height);
  if (ratio > 1){
    const elem = document.getElementById("style")
    elem.setAttribute("rel", "stylesheet");
  } else {
    const elem = document.getElementById("mobile_style");
    elem.setAttribute("rel", "stylesheet");
  }
  

})();
// This is a single line JS comment
/*
This is a comment that can span multiple lines 
- use comments to make your own notes!
*/
