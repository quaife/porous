function [xtar,ytar] = initialTracers(radii,centers,ntra,geom,...
    radiiBeans,centersBeans)
% find a set of intial tracer locations that are in the middle of the
% porous region
% keep xtar between 0.5 and 4.3
% keep ytar between 5 and 30

if nargin == 3
  geom = 'circles';
end


xmin = 2; xmax = 25; dx = xmax - xmin;
ymin = 0.1; ymax = 4.9; dy = ymax - ymin;
%xmin = 8.0; xmax = 8.30; dx = xmax - xmin;
%ymin = 2.1; ymax = 2.5; dy = ymax - ymin;
xtar = [];
ytar = [];
rng('shuffle');

if strcmp(geom(1:5),'circl')
  while numel(xtar) < ntra
    xtarPot = dx*rand(ceil(1.1*(ntra-numel(xtar))),1) + xmin;
    ytarPot = dy*rand(ceil(1.1*(ntra-numel(ytar))),1) + ymin;

    ncount = 0;
    nv = numel(radii);
    for k = 1:numel(xtarPot)
      if(~any((xtarPot(k) - centers(1:nv,1)).^2 + ...
            (ytarPot(k) - centers(1:nv,2)).^2 < radii(1:nv).^2))
        xtar = [xtar xtarPot(k)]; 
        ytar = [ytar ytarPot(k)]; 
      end
    end
  end
  xtar = xtar(1:ntra);
  ytar = ytar(1:ntra);
end

if strcmp(geom(1:5),'beans')
  while numel(xtar) < ntra
    indCutGrains = [360 139 403 195 116 253 323 398 263 364 362 ...
        223 306 349 431 179 1 396 109 315 418 277 295 374 310 183 ...
        260 289 447 185 313 338 143 404 307 128 280 347 324 406 420 ...
        316 287 416 434 370 436 312 432 402 415 384 369 405 408 333 ...
        417 410 373 351 443 371 353 279 321 298 411 422 392 388 231 ...
        424 345 445 452 236 429 217 400 401 439 379 234 380 309 412 ...
        461 361 302 433 464 221 346 435 459 414 419 391 423 383 442 ...
        219 428 437 375 465 305 451 457 441 462 395 440 367 463 430 ...
        444 448 446 427 449 460 409 454 450];

    indFullGrains = setdiff((1:465),[indCutGrains 421]);

    xtarPot = dx*rand(ceil(1.1*(ntra-numel(xtar))),1) + xmin;
    ytarPot = dy*rand(ceil(1.1*(ntra-numel(ytar))),1) + ymin;

    ncount = 0;
    nv = numel(radii);
    for k = 1:numel(xtarPot)
      if(~any((xtarPot(k) - centers(indFullGrains,1)).^2 + ...
            (ytarPot(k) - centers(indFullGrains,2)).^2 < ...
            radii(indFullGrains).^2) & ...
         (~any((xtarPot(k) - centers(indCutGrains,1)).^2 + ...
            (ytarPot(k) - centers(indCutGrains,2)).^2 < ...
            radii(indCutGrains).^2) | ...
          (any((xtarPot(k) - centersBeans(:,1)).^2 + ...
               (ytarPot(k) - centersBeans(:,2)).^2 < ...
               radiiBeans.^2)) ))
        xtar = [xtar xtarPot(k)]; 
        ytar = [ytar ytarPot(k)]; 
      end
    end
  end
  xtar = xtar(1:ntra);
  ytar = ytar(1:ntra);


%  xtar = [];
%  ytar = [];
%  xtar = dx*rand((ntra-numel(xtar)),1) + xmin;
%  ytar = dy*rand((ntra-numel(ytar)),1) + ymin;
%
end



