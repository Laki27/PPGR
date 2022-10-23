class Vector{
    constructor(x, y, z){
        this.x=x
        this.y=y
        this.z=z
    }
}

function crossProduct(v1, v2){
    return new Vector(v1.y * v2.z - v1.z * v2.y, 
                      v1.z * v2.x - v1.x * v2.z,
                      v1.x * v2.y - v1.y * v2.x)
}

function afinize(v){return new Vector(Math.floor(v.x/v.z), Math.floor(v.y/v.z), 1)}

const display_image=document.getElementById('display_image')
if(display_image==null){
    console.log("Nepostojece polje!!")
}
var width=0
var height=0
function podesi_visinu(){
    width=display_image.offsetWidth
    height=(9/16)*width
    display_image.style.height=height+'px'
    display_image.style.backgroundSize= `${width}px ${height}px` 
} 

const image_input=document.querySelector('#image_input')
var uploaded_image=""

var is_uploaded_image=false
i=0
image_input.addEventListener("change", function(){
    const reader = new FileReader()
    reader.addEventListener("load", ()=>{
        uploaded_image=reader.result
        document.querySelector("#display_image").style.backgroundImage= `url(${uploaded_image})`
    })
    reader.readAsDataURL(this.files[0]);
    const upload_icon = document.querySelector('img')
    upload_icon.style.display='none'
    is_uploaded_image=true
    i=0
})

setInterval(podesi_visinu, 1)

fields=document.querySelectorAll('.field')
vertex=document.querySelectorAll('.vertex')
result=document.getElementById('result')

points=[]

display_image.addEventListener('click', function(e){
        if(i<7 && is_uploaded_image){
            xOsa=e.clientX
            yOsa=e.clientY
            fields[i].value=`${xOsa}, ${yOsa}`
            fields[i].style.backgroundColor='blueviolet'
            points[i]=new Vector(xOsa, yOsa, 1)
            vertex[i].style.left=`${xOsa-8}px`
            vertex[i].style.top=`${yOsa-8}px`
            vertex[i].style.display='block'
            i++
            if(i==7){
                Xb=afinize(crossProduct(crossProduct(points[1], points[4]), crossProduct(points[0], points[3])))
                Yb=afinize(crossProduct(crossProduct(points[3], points[4]), crossProduct(points[5], points[6])))
                cilj=afinize(crossProduct(crossProduct(points[6], Xb), crossProduct(points[2], Yb)))
                res=document.querySelector('#res');
                res.innerHTML=`${cilj.x}, ${cilj.y}`
                result.style.left=`${cilj.x-9}px`
                result.style.top=`${cilj.y-9}px`
                result.style.display='block'
                i=0
            }
        }    
})