/** PGA3D.js
  * @author Adam thompson
  * @link   https://github.com/thehappyt/PGA3D.js
  */

function PGA3D() {
  const {abs,sqrt,pow,min,floor,cos,sin} = Math;

  var basis = ["1","e0","e1","e2","e3","e01","e02","e03","e12","e13","e23","e012","e013","e023","e123","e0123"];
  var grades = [0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,4];
  var grade_start = [0,1,5,11,15,16];
  var low = 0;
  
  var drm = [15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0];
  var drms = [1,-1,1,-1,1,1,-1,1,1,-1,1,-1,1,-1,1,1];
  var simplify = (s)=>{
    var sign=1, t=[]; s=s.match(/(\d)/g) || s;  if (!Array.isArray(s)) return s;
    for (var i=0; i<s.length-1; ) { if (s[i]>s[i+1]) { [s[i],s[i+1]] = [s[i+1],s[i]]; sign*=-1; i=0; } else ++i; }
    for (var i=0; i<s.length;) { if (s[i]===s[i+1]) { if (s[i]==0) sign=0; i+=2; } else t.push(s[i++]); }
    return ret=(sign==0)?'0':((sign==1)?'':'-')+(t.length?'e'+t.join(''):'1');
  };

  /// Generate multiplication tables for the outer and geometric products.
  var mulTable = basis.map(x=>basis.map(y=>(x==1)?y:(y==1)?x:simplify(x+y)));
  var gp=basis.map(x=>basis.map(y=>'0'));
  var cp=basis.map(x=>basis.map(y=>'0'));
  var cps=basis.map(x=>basis.map(y=>'0'));
  var op=basis.map(x=>basis.map(y=>'0'));
  var gpo={};
  basis.forEach((x,xi)=>basis.forEach((y,yi)=>{ var n = mulTable[xi][yi].replace(/^-/,''); if (!gpo[n]) gpo[n]=[]; gpo[n].push([xi,yi]); }));
  basis.forEach((o,oi)=>{
    gpo[o].forEach(([xi,yi])=>op[oi][xi]=(grades[oi]==grades[xi]+grades[yi])?((mulTable[xi][yi]=='0')?'0':((mulTable[xi][yi][0]!='-')?'':'-')+'b['+yi+']*this['+xi+']'):'0');
    gpo[o].forEach(([xi,yi])=>{
      gp[oi][xi] =((gp[oi][xi]=='0')?'':gp[oi][xi]+'+')   + ((mulTable[xi][yi]=='0')?'0':((mulTable[xi][yi][0]!='-')?'':'-')+'b['+yi+']*this['+xi+']');
      cp[oi][xi] =((cp[oi][xi]=='0')?'':cp[oi][xi]+'+')   + ((grades[oi]==grades[yi]-grades[xi])?gp[oi][xi]:'0');
      cps[oi][xi]=((cps[oi][xi]=='0')?'':cps[oi][xi]+'+') + ((grades[oi]==abs(grades[yi]-grades[xi]))?gp[oi][xi]:'0');
    });
  });

  var generator = class MultiVector extends (Float32Array) {
    constructor() { super(basis.length); return this; }
    
    Grade(grade,res) {
      res = res || new this.constructor();
      for (var i=0,l=res.length; i<l; i++) if (grades[i]==grade) res[i]=this[i]; else res[i]=0;
      return res;
    }
    
    Even(res) {
      res = res || new this.constructor();
      for (var i=0,l=res.length; i<l; i++) if (grades[i]%2==0) res[i]=this[i]; else res[i]=0;
      return res;
    }
    
    nVector(grade,...args) {
      this.set(args,grade_start[grade]);
      return this;
    }
    
    /// Fill in coordinates (accepts sequence of index,value as arguments)
    Coeff() {
      for (var i=0,l=arguments.length; i<l; i+=2) this[arguments[i]]=arguments[i+1];
      return this;
    }
    
    /// Negates specific grades (passed in as args)
    Map(res, ...a) {
      for (var i=0, l=res.length; i<l; i++) res[i] = (~a.indexOf(grades[i]))?-this[i]:this[i];
      return res;
    }
    
    /// Returns the vector grade only.
    get Vector () {
      return this.slice(grade_start[1],grade_start[2]);
    };
    
    toString() {
      var res=[];
      for (var i=0; i<basis.length; i++) if (Math.abs(this[i])>1e-10) res.push(((this[i]==1)&&i?'':((this[i]==-1)&&i)?'-':(this[i].toFixed(10)*1))+(i==0?'':basis[i].replace('e','e_')));
      return res.join('+').replace(/\+-/g,'-')||'0';
    }
    
    /// Reversion, Involutions, Conjugation for any number of grades, component acces shortcuts.
    get Negative (){
      var res = new this.constructor();
      for (var i=0; i<this.length; i++) res[i] = -this[i];
      return res;
    };
    get Reverse (){
      var res = new this.constructor();
      for (var i=0; i<this.length; i++) res[i]= this[i]*[1,1,-1,-1][grades[i]%4];
      return res;
    };
    get Involute (){
      var res = new this.constructor();
      for (var i=0; i<this.length; i++) res[i]= this[i]*[1,-1,1,-1][grades[i]%4];
      return res;
    };
    get Conjugate (){
      var res = new this.constructor();
      for (var i=0; i<this.length; i++) res[i]= this[i]*[1,-1,-1,1][grades[i]%4];
      return res;
    };

    /// The Dual, Length, non-metric length and normalized getters.
    get Dual (){ return this.map((x,i,a)=>a[drm[i]]*drms[i]); };
    get UnDual (){ return this.map((x,i,a)=>a[drm[i]]*drms[a.length-i-1]); };
    get Length (){ return Math.sqrt(Math.abs(this.Mul(this.Conjugate).s)); };
    get VLength (){
      var res = 0;
      for (var i=0; i<this.length; i++) res += this[i]*this[i];
      return Math.sqrt(res);
    };
    get Normalized (){ 
      var sq = this.Mul(this.Reverse), [s,t] = [sq[0], sq[15]];
      if (s==0) return this;
      sq[0] = 1/Math.sqrt(s); sq[15] = -t/(2*Math.pow(Math.sqrt(s),3));
      return this.Mul(sq);
    };
  }
  
  /// Convert symbolic matrices to code. (skipping zero's on dot and wedge matrices).
  /// These all do straightforward string fiddling. If the 'mix' option is set they reference basis components using e.g. '.e1' instead of eg '[3]' .. so that
  /// it will work for elements of subalgebras etc.
  generator.prototype.Add   = new Function('b,res','res=res||new this.constructor();\n'+basis.map((x,xi)=>'res['+xi+']=b['+xi+']+this['+xi+']').join(';\n')+';\nreturn res')
  generator.prototype.Scale = new Function('b,res','res=res||new this.constructor();\n'+basis.map((x,xi)=>'res['+xi+']=b*this['+xi+']').join(';\n')+';\nreturn res')
  generator.prototype.Sub   = new Function('b,res','res=res||new this.constructor();\n'+basis.map((x,xi)=>'res['+xi+']=this['+xi+']-b['+xi+']').join(';\n')+';\nreturn res')
  //generator.prototype.Mul   = new Function('b,res','res=res||new this.constructor();\n'+gp.map((r,ri)=>'res['+ri+']='+r.join('+').replace(/\+\-/g,'-').replace(/(\w*?)\[(.*?)\]/g,(a,b,c)=>options.mix?'('+b+'.'+(c|0?basis[c]:'s')+'||0)':a).replace(/\+0/g,'')+';').join('\n')+'\nreturn res;');
  //generator.prototype.LDot  = new Function('b,res','res=res||new this.constructor();\n'+cp.map((r,ri)=>'res['+ri+']='+r.join('+').replace(/\+\-/g,'-').replace(/\+0/g,'').replace(/(\w*?)\[(.*?)\]/g,(a,b,c)=>options.mix?'('+b+'.'+(c|0?basis[c]:'s')+'||0)':a)+';').join('\n')+'\nreturn res;');
  //generator.prototype.Dot   = new Function('b,res','res=res||new this.constructor();\n'+cps.map((r,ri)=>'res['+ri+']='+r.join('+').replace(/\+\-/g,'-').replace(/\+0/g,'').replace(/(\w*?)\[(.*?)\]/g,(a,b,c)=>options.mix?'('+b+'.'+(c|0?basis[c]:'s')+'||0)':a)+';').join('\n')+'\nreturn res;');
  //generator.prototype.Wedge = new Function('b,res','res=res||new this.constructor();\n'+op.map((r,ri)=>'res['+ri+']='+r.join('+').replace(/\+\-/g,'-').replace(/\+0/g,'').replace(/(\w*?)\[(.*?)\]/g,(a,b,c)=>options.mix?'('+b+'.'+(c|0?basis[c]:'s')+'||0)':a)+';').join('\n')+'\nreturn res;');
  // generator.prototype.Vee   = new Function('b,res','res=res||new this.constructor();\n'+op.map((r,ri)=>'res['+drm[ri]+']='+r.map(x=>x.replace(/\[(.*?)\]/g,function(a,b){return '['+(drm[b|0])+']'})).join('+').replace(/\+\-/g,'-').replace(/\+0/g,'').replace(/(\w*?)\[(.*?)\]/g,(a,b,c)=>options.mix?'('+b+'.'+(c|0?basis[c]:'s')+'||0)':a)+';').join('\n')+'\nreturn res;');
  /// Conforms to the new Chapter 11 now.
  //generator.prototype.Vee   = new Function('b,res',('res=res||new this.constructor();\n'+op.map((r,ri)=>'res['+drm[ri]+']='+drms[ri]+'*('+r.map(x=>x.replace(/\[(.*?)\]/g,function(a,b){return '['+(drm[b|0])+']'+(drms[b|0]>0?"":"*-1")})).join('+').replace(/\+\-/g,'-').replace(/\+0/g,'').replace(/(\w*?)\[(.*?)\]/g,(a,b,c)=>options.mix?'('+b+'.'+(c|0?basis[c]:'s')+'||0)':a)+');').join('\n')+'\nreturn res;').replace(/(b\[)|(this\[)/g,a=>a=='b['?'this[':'b['));
  //generator.prototype.eigenValues = eigenValues;
  
  /// Add getter and setters for the basis vectors/bivectors etc ..
  basis.forEach((b,i)=>Object.defineProperty(generator.prototype, i?b:'s', {
    configurable: true, get(){ return this[i] }, set(x){ this[i]=x; }
  }));
  
}


/*
    var simplify_bits = (A,B,p2)=>{
        var n=p2||(p+q+r), t=0, ab=A&B, res=A^B;
        if (ab&((1<<r)-1)) return [0,0];
        while (n--) t^=(A=A>>1);
        t&=B;
        t^=ab>>(p+r);
        t^=t>>16;
        t^=t>>8;
        t^=t>>4;
        return [1-2*(27030>>(t&15)&1),res];
      },
      bc = (v)=>{
        v=v-((v>>1)& 0x55555555);
        v=(v&0x33333333)+((v>>2)&0x33333333);
        var c=((v+(v>>4)&0xF0F0F0F)*0x1010101)>>24;
        return c
      };
*/
