// sourcery.js — static JavaScript for all Sourcery pages.
// Object-specific values are injected by the server into window.sourceryConfig
// before this file is loaded (on the source page).

// ── Index page: collapsible fieldset legends ────────────────────────────────

$(function () {
    $("legend:has(span.toggle)").css("cursor", "pointer").click(function () {
        var legend = $(this);
        var value = $(this).children("span.toggle").html();
        value = (value === "hide") ? "expand" : "hide";
        $(this).siblings().slideToggle(150, function () {
            legend.children("span.toggle").html(value);
        });
    });
});

$(document).ready(function () {
    $("legend span.toggle").each(function () {
        if ($(this).html() === "expand") {
            $(this).html("hide");
            $(this).closest("legend").click();
        }
    });
});

// ── Source page: image controls ─────────────────────────────────────────────
// Only runs when sourceryConfig is present (i.e., on the source page).

if (typeof sourceryConfig !== "undefined") {

    function printValue(sliderID, textbox) {
        document.getElementById(textbox).value = document.getElementById(sliderID).value;
    }

    function updateSlider(value) {
        var slider = document.getElementById("sizeSlider");
        slider.value = value;
        printValue("sizeSlider", "sizeSliderValue");
    }

    function parseImageTypeValue(imageType, index) {
        return imageType.split(";")[index];
    }

    window.onload = function () {
        printValue("sizeSlider", "sizeSliderValue");
        printValue("gammaSlider", "gammaSliderValue");
    };

    // Load the initial image on page load.
    $(document).ready(function () {
        $.post("makePlotFromJPEG",
            {
                name:           sourceryConfig.objectName,
                RADeg:          sourceryConfig.RADeg,
                decDeg:         sourceryConfig.decDeg,
                surveyLabel:    parseImageTypeValue($("input:radio[name=imageType]:checked").val(), 0),
                plotNEDObjects: $("input:checkbox[name=plotNEDObjects]").prop("checked"),
                plotSpecObjects:$("input:checkbox[name=plotSpecObjects]").prop("checked"),
                plotSourcePos:  $("input:checkbox[name=plotSourcePos]").prop("checked"),
                plotXMatch:     $("input:checkbox[name=plotXMatch]").prop("checked"),
                plotContours:   $("input:checkbox[name=plotContours]").prop("checked"),
                showAxes:       $("input:checkbox[name=showAxes]").prop("checked"),
                clipSizeArcmin: sourceryConfig.defaultClipSizeArcmin,
                gamma:          $("#gamma").val(),
                plotRedshift:   $("input:checkbox[name=plotRedshift]").prop("checked"),
                redshift:       sourceryConfig.redshift
            },
            function (data) {
                $("#imagePlot").html(
                    '<img src="data:image/jpg;base64,' + data +
                    '" align="middle" border=0 width="' + sourceryConfig.plotDisplayWidthPix + '"/>'
                );
            }
        );
    });

    // Resubmit the form when any checkbox or radio button changes.
    $(document).ready(function () {
        $(':checkbox').change(function () { $("#imageForm").submit(); });
        $('input:radio[name=imageType]').change(function () { $("#imageForm").submit(); });
    });

    // Handle form submission (image reload).
    $(function () {
        $("#imageForm").submit(function () {
            var maxSizeArcmin = parseImageTypeValue($("input:radio[name=imageType]:checked").val(), 1);
            if (maxSizeArcmin > sourceryConfig.defaultClipSizeArcmin) {
                updateSlider(parseImageTypeValue($("input:radio[name=imageType]:checked").val(), 1));
            }
            if (parseFloat($("#sizeSliderValue").val()) > maxSizeArcmin) {
                updateSlider(maxSizeArcmin);
            }
            $.post("makePlotFromJPEG",
                {
                    name:           sourceryConfig.objectName,
                    RADeg:          sourceryConfig.RADeg,
                    decDeg:         sourceryConfig.decDeg,
                    surveyLabel:    parseImageTypeValue($("input:radio[name=imageType]:checked").val(), 0),
                    plotNEDObjects: $("input:checkbox[name=plotNEDObjects]").prop("checked"),
                    plotSpecObjects:$("input:checkbox[name=plotSpecObjects]").prop("checked"),
                    plotSourcePos:  $("input:checkbox[name=plotSourcePos]").prop("checked"),
                    plotXMatch:     $("input:checkbox[name=plotXMatch]").prop("checked"),
                    plotContours:   $("input:checkbox[name=plotContours]").prop("checked"),
                    showAxes:       $("input:checkbox[name=showAxes]").prop("checked"),
                    clipSizeArcmin: $("#sizeSliderValue").val(),
                    gamma:          $("#gammaSliderValue").val(),
                    plotRedshift:   $("input:checkbox[name=plotRedshift]").prop("checked"),
                    redshift:       sourceryConfig.redshift
                },
                function (data) {
                    $("#imagePlot").html(
                        '<img src="data:image/jpg;base64,' + data +
                        '" align="middle" border=0 width="' + sourceryConfig.plotDisplayWidthPix + '"/>'
                    );
                }
            );
            return false;
        });
    });

    // Clickable-coordinate handler — only registered when fields are configured.
    if (sourceryConfig.clickableRAField) {
        $(document).ready(function () {
            $("#imagePlot").on("click", function (event) {
                var x = event.pageX - this.offsetLeft;
                var y = event.pageY - this.offsetTop;
                var width  = $("#imagePlot").width();
                var height = $("#imagePlot").height();
                var CDELT1 = ($("#sizeSliderValue").val() / 60.0) / width;
                var CDELT2 = ($("#sizeSliderValue").val() / 60.0) / height;
                var CRPIX1 = (width  / 2.0) + 1;
                var CRPIX2 = (height / 2.0) + 1;
                var CRVAL1 = sourceryConfig.RADeg;
                var CRVAL2 = sourceryConfig.decDeg;

                x = width  - x;
                y = height - y;

                var deg2Rad = Math.PI / 180;
                var xint = deg2Rad * (CDELT1 * (x - CRPIX1));
                var yint = deg2Rad * (CDELT2 * (y - CRPIX2));

                var Rtheta  = Math.sqrt(Math.pow(xint, 2) + Math.pow(yint, 2));
                var phi     = Math.atan2(xint, -yint);
                var theta   = Math.atan(1.0 / Rtheta);

                var alpha_p = deg2Rad * CRVAL1;
                var delta_p = deg2Rad * CRVAL2;
                var phi_p   = deg2Rad * 180.0;

                var ra  = alpha_p + Math.atan2(
                    -Math.cos(theta) * Math.sin(phi - phi_p),
                    Math.sin(theta) * Math.cos(delta_p) - Math.cos(theta) * Math.sin(delta_p) * Math.cos(phi - phi_p)
                );
                var dec = Math.asin(
                    Math.sin(theta) * Math.sin(delta_p) +
                    Math.cos(theta) * Math.cos(delta_p) * Math.cos(phi - phi_p)
                );

                ra  = (180 / Math.PI) * ra;
                dec = (180 / Math.PI) * dec;

                if ($("input:checkbox[name=showAxes]").prop("checked")) {
                    alert("Coords not recorded when showing coordinate axes — deselect 'Show coordinate axes'");
                } else {
                    $("[name=" + sourceryConfig.clickableRAField  + "]").val(ra);
                    $("[name=" + sourceryConfig.clickableDecField + "]").val(dec);
                }
            });
        });
    }

} // end if (typeof sourceryConfig !== "undefined")
